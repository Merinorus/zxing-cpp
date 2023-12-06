/*
 * Copyright 2023 ZXing authors
 * Copyright 2023 Antoine Mérino
 */
// SPDX-License-Identifier: Apache-2.0

#include "DecodeHints.h"
#include "GTIN.h"
#include "ODDXFilmEdgeReader.h"
#include "Result.h"
#include "ZXAlgorithms.h"

#include <array>
#include <iostream>
#include <bitset>
#include <iomanip>
#include<numeric>
namespace ZXing::OneD {

// Detection is made from center to bottom.
// We ensure the clock signal is decoded before the data signal to avoid false positives.
// They are two version of a DX Edge codes : without half-frame information and with half-frame information.
// The clock signal is longer if the DX code contains the half-frame information (more recent version)
const int DATA_START_PATTERN_SIZE = 5;
constexpr auto CLOCK_PATTERN_HF = FixedPattern<25, 31>   {5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3};
constexpr auto CLOCK_PATTERN_NO_HF = FixedPattern<17, 23>{5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3};
constexpr auto DATA_START_PATTERN_ = FixedPattern<5, 5> {1, 1, 1, 1, 1};
constexpr auto DATA_STOP_PATTERN_ = FixedPattern<3, 3> {1, 1, 1};

/**
 * @brief Parse a part of a vector of bits (boolean) to a decimal number.
 * Eg: {1, 1, 0} -> 6.
 * @param begin begin of the vector's part to be parsed
 * @param end end of the vector's part to be parsed
 * @return The decimal value of the parsed part
 */
int binaryToDecimal(std::vector<bool>::iterator begin, std::vector<bool>::iterator end)
{
	int retval = 0;
	auto i = std::distance(begin, end) - 1;
	for (std::vector<bool>::iterator it = begin; it != end; it++, i--) {
		retval += (*it ? (1 << i) : 0);
	}
	return retval;
}

// To avoid many false positives,
// the clock signal must be found to attempt to decode a data signal.
// We ensure the data signal starts below a clock.
// We accept a tolerance margin,
// ie. the signal may start a few pixels before or after the clock on the X-axis.
struct DXFEState : public RowReader::DecodingState
{
	bool foundClock = false;
	bool containsHFNumber = false; // Clock signal (thus data signal) with half-frame number (longer version)
	int clockXStart = 0; // Beginning of the clock signal on the X-axis, in pixels
	int clockXStop = 0; // End of the clock signal on the X-axis, in pixels 
	int pixelTolerance = 0; // Pixel tolerance will be set depending of the length of the clock signal (in pixels)
};

Result DXFilmEdgeReader::decodePattern(int rowNumber, PatternView& next, std::unique_ptr<DecodingState>& state) const
{
	// Retrieve the decoding state to check if a clock signal has already been found before.
	// We check also if it contains the half-frame number, since it affects the signal structure.
	if (!state)
		state.reset(new DXFEState);
	auto& foundClock = static_cast<DXFEState*>(state.get())->foundClock;
	auto& clockXStart = static_cast<DXFEState*>(state.get())->clockXStart;
	auto& clockXStop = static_cast<DXFEState*>(state.get())->clockXStop;
	auto& pixelTolerance = static_cast<DXFEState*>(state.get())->pixelTolerance;
	auto& containsHFNumber = static_cast<DXFEState*>(state.get())->containsHFNumber;

	// TODO: what is this measure exactly?
	const int minCharCount = 5;

	// Minimum "allowed "white" zone to the left and the right sides of the signal.
	// The margin is low because of potential sprocket holes close to the signal.
	const float minQuietZone = 0.4;

	// Adjust the pixel shift tolerance between the data signal and the clock signal.
	// 1 means the signal can be shifted up to one bar to the left or the right.
	const float pixelToleranceRatio = 0.4;

	// Since we parse the image from the center to the edges (top & bottom)
	// The clock should be found before the data
	// So we always make sure we have found the clock before parsing the data signal.
	if (!foundClock) {
		auto clock = FindLeftGuard(next, minCharCount, CLOCK_PATTERN_HF, minQuietZone);
		if (clock.isValid()) {
			foundClock = true;
			containsHFNumber = true;
		} else {
			clock = FindLeftGuard(next, minCharCount, CLOCK_PATTERN_NO_HF, minQuietZone);
			if (clock.isValid())
				foundClock = true;
		}
		if (foundClock) {
			clockXStart = clock.pixelsInFront();
			clockXStop = clock.pixelsTillEnd();
			pixelTolerance = (clock.pixelsTillEnd() - clock.pixelsInFront()) / (containsHFNumber ? 31 : 23) * pixelToleranceRatio;
		}
	}

	if (!foundClock)
		// Clock not found. Aborting the decoding of the current row.
		return {};



	// Now that we found the clock, attempt to decode the data signal.
	// Start by finding the data start pattern.
	next = FindLeftGuard(next, minCharCount, DATA_START_PATTERN_, minQuietZone);
	if (!next.isValid())
		return {};

	int xStart = next.pixelsInFront();

	// The found data signal must be below the clock signal, otherwise we abort the decoding (potential false positive)
	if (!(xStart - pixelTolerance < clockXStart && xStart + pixelTolerance > clockXStart)) {
		return {};
	}

	//  Compute the length of a bar
	// It may be greater than 1 depending on what have been found in the raw signal
	auto per_bar_raw_width = *next.data();

	// Skip the data start pattern (black white black white black)
	// The first signal bar is always white
	// (separation between the start pattern and the product number) 
	next.shift(DATA_START_PATTERN_SIZE);

	//auto signal_begin = next;
	int signal_length = 0;

	if (!next.isValid())
		return {};

	std::vector<bool> signal_data;
	// They are two possible data signal lengths (with or without half-frame information)
	signal_data.reserve(containsHFNumber ? 23 : 15);
	bool current_signal_is_black = false;

	// Populate a vector of booleans to represent the bits. true = black, false = white.
	// We start the parsing just after the data start signal.
	// The first bit is always a white bar (we include the separator just after the start pattern)
	// Eg: {3, 1, 2} -> {0, 0, 0, 1, 0, 0}
	while (signal_length < (containsHFNumber ? 23 : 15)) {
		if (!next.isValid())
			return {};

		if (*next.data() == 0) {
			// Zero means we can't conclude on black or white bar. Abort the decoding.
			return {};
		}

		auto current_bar_width =
			*next.data() / per_bar_raw_width + (*next.data() % per_bar_raw_width >= (per_bar_raw_width / 2) ? 1 : 0);

		signal_length += current_bar_width;
		while (current_bar_width > 0 && signal_data.size() < (containsHFNumber ? 23 : 15)) {
			signal_data.push_back(current_signal_is_black);
			--current_bar_width;
		}
		current_signal_is_black = !current_signal_is_black;
		next.shift(1);
	}

	// Check the signal length
	if (signal_length != (containsHFNumber ? 23 : 15))
		return {};
	
	// Check there is the Stop pattern at the end of the data signal
	next = next.subView(0, 3);
	if (!IsRightGuard(next, DATA_STOP_PATTERN_, minQuietZone)) {
		return {};
	}

	// Check the data signal has been fully parsed
	if (containsHFNumber && signal_data.size() < 23)
		return {};

	if (!containsHFNumber && signal_data.size() < 15)
		return {};

	// The following bits are always white (=false), they are separators.
	if (signal_data.at(0) || signal_data.at(8))
		return {};

	if (containsHFNumber && (signal_data.at(20) || signal_data.at(22)))
		return {};

	if (!containsHFNumber && (signal_data.at(8) || signal_data.at(14)))
		return {};

	// Check we didn't just parse the clock signal instead of the data signal.
	// It might be parsed accidentally when the minQuietZone is too small.
	if (minQuietZone <= 1 && containsHFNumber) {
		std::vector<bool> signal_clock = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0};
		if (signal_clock == signal_data)
			return {};
	}

	// Check the parity bit
	auto signal_sum = std::accumulate(signal_data.begin(), signal_data.end()-2, 0);
	auto parity_bit = *(signal_data.end() - 2);
	if (signal_sum % 2 != parity_bit)
		return {};

	// Compute the DX 1 number (product number)
	auto product_number = binaryToDecimal(signal_data.begin() + 1, signal_data.begin() + 8);
	if (!product_number)
		return {};

	// Compute the DX 2 number (generation number)
	auto generation_number = binaryToDecimal(signal_data.begin() + 9, signal_data.begin() + 13);

	// Generate the textual representation.
	// Eg: 115-10/11A means: DX1 = 115, DX2 = 10, Frame number = 11A
	std::string txt;
	txt.reserve(10);
	txt = std::to_string(product_number) + "-" + std::to_string(generation_number);
	if (containsHFNumber) {
		auto half_frame_number = binaryToDecimal(signal_data.begin() + 13, signal_data.begin() + 20);
		txt += "/" + std::to_string(half_frame_number / 2);
		if (half_frame_number % 2)
			txt += "A";
	}

	Error error;

	// TODO is it required?
	// AFAIK The DX Edge barcode doesn't follow any symbology identifier.
	SymbologyIdentifier symbologyIdentifier = {'I', '0'}; // No check character validation ?

	int xStop = next.pixelsTillEnd();

	// The found data signal must be below the clock signal, otherwise we abort the decoding (potential false positive)
	if (!(xStop - pixelTolerance < clockXStop && xStop + pixelTolerance > clockXStop)) {
		return {};
	}

	return Result(txt, rowNumber, xStart, xStop, BarcodeFormat::DXFilmEdge, symbologyIdentifier, error);
}

} // namespace ZXing::OneD
