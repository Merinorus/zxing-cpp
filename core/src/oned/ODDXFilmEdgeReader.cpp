/*
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
#include <numeric>
#include <set>
namespace ZXing::OneD {

// Detection is made from center to bottom.
// We ensure the clock signal is decoded before the data signal to avoid false positives.
// They are two version of a DX Edge codes : without half-frame information and with half-frame information.
// The clock signal is longer if the DX code contains the half-frame information (more recent version)
const int CLOCK_PATTERN_LENGTH_HF = 31;
const int CLOCK_PATTERN_LENGTH_NO_HF = 23;
const int DATA_START_PATTERN_SIZE = 5;
constexpr auto CLOCK_PATTERN_COMMON = FixedPattern<15, 20> {5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr auto CLOCK_PATTERN_HF =
	FixedPattern<25, CLOCK_PATTERN_LENGTH_HF>{5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3};
constexpr auto CLOCK_PATTERN_NO_HF = FixedPattern<17, CLOCK_PATTERN_LENGTH_NO_HF>{5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3};
constexpr auto DATA_START_PATTERN_ = FixedPattern<5, 5>    {1, 1, 1, 1, 1};
constexpr auto DATA_STOP_PATTERN_ = FixedPattern<3, 3>     {1, 1, 1};

// Signal data length, without the start and stop patterns
const int DATA_LENGTH_HF = 23;
const int DATA_LENGTH_NO_HF = 15;

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

class Clock
{
public:
	int rowNumber = 0;
	bool containsHFNumber = false; // Clock signal (thus data signal) with half-frame number (longer version)
	int xStart = 0; // Beginning of the clock signal on the X-axis, in pixels
	int xStop = 0; // End of the clock signal on the X-axis, in pixels
	int pixelTolerance = 0; // Pixel tolerance will be set depending of the length of the clock signal (in pixels)

	bool operator<(const int x) const
	{
		return xStart < x;
	}

	bool operator < (const Clock& other) const {
		return xStart < other.xStart;
	}

	// We assume two clock are the same when they start in about the same X position,
	// even if they are different clocks (stop at different position or different type)
	// Only the more recent clock is kept.
	bool xStartInRange(const Clock& other) const {
		auto tolerance = std::max(pixelTolerance, other.pixelTolerance);
		return (xStart - tolerance) <= other.xStart && (xStart + tolerance) >= other.xStart;
	}

	bool xStartInRange(const int x) const
	{ return (xStart - pixelTolerance) <= x && (xStart + pixelTolerance) >= x;
	}

	bool xStopInRange(const int x) const { return (xStop - pixelTolerance) <= x && (xStop + pixelTolerance) >= x; }

};

class ClockSet : public std::set<Clock, std::less<>>
{
public:

	// Return the closest clock on the X axis
	const ClockSet::iterator closest_element(const int x) {
		const ClockSet::iterator it = lower_bound(x);
		if (it == begin())
			return it;

		const ClockSet::iterator prev_it = std::prev(it);
		return (it == end() || x - prev_it->xStart <= it->xStart - x) ? prev_it : it;
	}
};


// To avoid many false positives,
// the clock signal must be found to attempt to decode a data signal.
// We ensure the data signal starts below a clock.
// We accept a tolerance margin,
// ie. the signal may start a few pixels before or after the clock on the X-axis.
struct DXFEState : public RowReader::DecodingState
{
	ClockSet allClocks;
};

template <typename Set>
//Helper to find the closest element in a set 
auto closest_element(Set& set, const typename Set::value_type& value) -> decltype(set.begin())
{
	const auto it = set.lower_bound(value);
	if (it == set.begin())
		return it;

	const auto prev_it = std::prev(it);
	return (it == set.end() || value - *prev_it <= *it - value) ? prev_it : it;
}


void findClock(int rowNumber, PatternView& next, ClockSet& allClocks)
{
	// Minimum "allowed "white" zone to the left and the right sides of the clock signal.
	const float minClockNoHFQuietZone = 2;
	const float minClockHFQuietZone = 1;

	// Adjust the pixel shift tolerance between the data signal and the clock signal.
	// 1 means the signal can be shifted up to one bar to the left or the right.
	const float pixelToleranceRatio = 0.5;

	// Before detecting any clock,
	// try to detect the common pattern between all types of clocks.
	// This avoid doing two detections at each interations instead of one,
	// when they is no DX Edge code to detect.
	auto common_clock_pattern =
		FindLeftGuard(next, CLOCK_PATTERN_COMMON.size(), CLOCK_PATTERN_COMMON, std::min(minClockNoHFQuietZone, minClockHFQuietZone));
	if (common_clock_pattern.isValid()) {
		bool foundClock = false;
		bool containsHFNumber = false;
		auto clock_pattern = FindLeftGuard(next, CLOCK_PATTERN_HF.size(), CLOCK_PATTERN_HF, minClockHFQuietZone);
		if (clock_pattern.isValid()) {
			foundClock = true;
			containsHFNumber = true;
		} else {
			clock_pattern = FindLeftGuard(next, CLOCK_PATTERN_NO_HF.size(), CLOCK_PATTERN_NO_HF, minClockNoHFQuietZone);
			if (clock_pattern.isValid())
				foundClock = true;
		}
		if (foundClock) {
			Clock clock;
			clock.rowNumber = rowNumber;
			clock.containsHFNumber = containsHFNumber;
			clock.xStart = clock_pattern.pixelsInFront();
			clock.xStop = clock_pattern.pixelsTillEnd();
			clock.pixelTolerance = (clock_pattern.pixelsTillEnd() - clock_pattern.pixelsInFront())
								   / (containsHFNumber ? CLOCK_PATTERN_LENGTH_HF : CLOCK_PATTERN_LENGTH_NO_HF) * pixelToleranceRatio;
			// Check if the clock was not already found
			auto closestClock = allClocks.closest_element(clock.xStart);
			if (closestClock != allClocks.end() && clock.xStartInRange(*closestClock)) {
				// If the clock was already found, replace with the new coordinates
				// This allow to find the signal data more precisely if the image is skewed
				allClocks.erase(closestClock);
				allClocks.insert(clock);
			} else {
				allClocks.insert(clock);
			}
		}
	}



}

Result DXFilmEdgeReader::decodePattern(int rowNumber, PatternView& next, std::unique_ptr<DecodingState>& state) const
{
	// Retrieve the decoding state to check if a clock signal has already been found before.
	// We check also if it contains the half-frame number, since it affects the signal structure.
	if (!state)
		state.reset(new DXFEState);
	auto& allClocks = static_cast<DXFEState*>(state.get())->allClocks;



	// Minimum "allowed "white" zone to the left and the right sides of the data signal.
	// We allow a smaller quiet zone, ie improve detection at risk of getting false positives,
	// because the risk is greatly reduced when we check we found the clock before the signal.
	const float minDataQuietZone = 0.5;

	// TODO Separate in two separate functions? lack of separation of concerns
	findClock(rowNumber, next, allClocks);

	// We should find at least one clock before attempting to decode the data signal.
	if (allClocks.empty())
		return {};


	// Now that we found at least a clock, attempt to decode the data signal.
	// Start by finding the data start pattern.
	next = FindLeftGuard(next, DATA_START_PATTERN_.size(), DATA_START_PATTERN_, minDataQuietZone);
	if (!next.isValid())
		return {};

	int xStart = next.pixelsInFront();

	// The found data signal must be below the clock signal, otherwise we abort the decoding (potential false positive)
	auto closestClock = allClocks.closest_element(xStart);
	if (!closestClock->xStartInRange(xStart)) {
		return {};
	}

	// Avoid decoding a signal found at the top of the clock
	// (might happen when stacking two films one of top of the other)
	if (closestClock->rowNumber > rowNumber) {
		return {};
	}

	//  Compute the length of a bar
	// It may be greater than 1 depending on what have been found in the raw signal
	auto per_bar_raw_width = *next.data();

	// Skip the data start pattern (black, white, black, white, black)
	// The first signal bar is always white: this is the
	// separation between the start pattern and the product number)
	next.shift(DATA_START_PATTERN_SIZE);

	if (!next.isValid())
		return {};

	std::vector<bool> signal_data;

	// They are two possible data signal lengths (with or without half-frame information)
	signal_data.reserve(closestClock->containsHFNumber ? DATA_LENGTH_HF : DATA_LENGTH_NO_HF);

	// Populate a vector of booleans to represent the bits. true = black, false = white.
	// We start the parsing just after the data start signal.
	// The first bit is always a white bar (we include the separator just after the start pattern)
	// Eg: {3, 1, 2} -> {0, 0, 0, 1, 0, 0}
	int signal_length = 0;
	bool current_bar_is_black = false; // the current bar is white
	while (signal_length < (closestClock->containsHFNumber ? DATA_LENGTH_HF : DATA_LENGTH_NO_HF)) {
		if (!next.isValid())
			return {};

		// Zero means we can't conclude on black or white bar. Abort the decoding.
		if (*next.data() == 0) {
			return {};
		}

		// Adjust the current bar according to the computed ratio above.
		// When the raw result is not exact (between two bars),
		// we round the bar size to the nearest integer.
		auto current_bar_width =
			*next.data() / per_bar_raw_width + (*next.data() % per_bar_raw_width >= (per_bar_raw_width / 2) ? 1 : 0);

		signal_length += current_bar_width;

		// Extend the bit array according to the current bar length.
		// Eg: one white bars -> {0}, three black bars -> {1, 1, 1}
		while (current_bar_width > 0
			   && signal_data.size() < (closestClock->containsHFNumber ? DATA_LENGTH_HF : DATA_LENGTH_NO_HF)) {
			signal_data.push_back(current_bar_is_black);
			--current_bar_width;
		}

		// Iterate to the next data signal bar (the color is inverted)
		current_bar_is_black = !current_bar_is_black;
		next.shift(1);
	}


	// Check the signal length
	if (signal_length != (closestClock->containsHFNumber ? DATA_LENGTH_HF : DATA_LENGTH_NO_HF))
		return {};
	
	// Check there is the Stop pattern at the end of the data signal
	next = next.subView(0, 3);
	if (!IsRightGuard(next, DATA_STOP_PATTERN_, minDataQuietZone)) {
		return {};
	}
	// Check the data signal has been fully parsed
	if (closestClock->containsHFNumber && signal_data.size() < DATA_LENGTH_HF)
		return {};
	if (!closestClock->containsHFNumber && signal_data.size() < DATA_LENGTH_NO_HF)
		return {};

	// The following bits are always white (=false), they are separators.
	if (signal_data.at(0) || signal_data.at(8))
		return {};
	if (closestClock->containsHFNumber && (signal_data.at(20) || signal_data.at(22)))
		return {};
	if (!closestClock->containsHFNumber && (signal_data.at(8) || signal_data.at(14)))
		return {};

	// Check we didn't just parse the clock signal instead of the data signal.
	// It might be parsed accidentally when the minQuietZone is too small.
	if (minDataQuietZone <= 1 && closestClock->containsHFNumber) {
		std::vector<bool> signal_clock = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0};
		if (signal_clock == signal_data)
			return {};
	}

	// Check the parity bit
	auto signal_sum = std::accumulate(signal_data.begin(), signal_data.end()-2, 0);
	auto parity_bit = *(signal_data.end() - 2);
	if (signal_sum % 2 != (int)parity_bit)
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
	if (closestClock->containsHFNumber) {
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
	if (!closestClock->xStopInRange(xStop)) {
		return {};
	}

	// Update the clock X coordinates with the latest corresponding data signal
	// This may improve signal detection for next row iterations
	if (closestClock->xStop != xStop || closestClock->xStart != xStart) {
		Clock clock(*closestClock);
		clock.xStart = xStart;
		clock.xStop = xStop;
		allClocks.erase(closestClock);
		allClocks.insert(clock);
	}

	return Result(txt, rowNumber, xStart, xStop, BarcodeFormat::DXFilmEdge, symbologyIdentifier, error);
}

} // namespace ZXing::OneD
