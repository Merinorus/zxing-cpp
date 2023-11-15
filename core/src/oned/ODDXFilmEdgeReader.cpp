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

// Pure binary translation. Nothing to translate actually. Should be simplified.
//static const char ALPHABET[] = "01";
//static const int CHARACTER_ENCODINGS[] = {0b0, 0b1};


constexpr auto START_PATTERN_ = FixedPattern<5, 5>{1, 1, 1, 1, 1};
constexpr auto STOP_PATTERN_ = FixedPattern<3, 3>{1, 1, 1};
//constexpr auto START_PATTERN_ = FixedPattern<5, 10>{2, 2, 2, 2, 2};
//constexpr auto STOP_PATTERN_ = FixedPattern<3, 6>{2, 2, 2};

int fromVector(std::vector<bool>::iterator begin, std::vector<bool>::iterator end)
{
	//std::cout << "fromVector: ";
	int retval = 0;
	auto i = std::distance(begin, end) - 1;
	for (std::vector<bool>::iterator it = begin; it != end; it++, i--) {

		//std::cout << "i=" << i << std::endl;
		//std::cout << "id=" << (int)*it << std::endl;
		//std::cout << "-> increment=" << (*it ? (1 << i) : 0) << std::endl;
		retval += (*it ? (1 << i) : 0);
		//if (*it) {
		//	retval += 1 << i;
		//}
		/*retval |= 1 << (int)*it;*/
	}
	//std::cout << std::endl;
	return retval;
}

Result DXFilmEdgeReader::decodePattern(int rowNumber, PatternView& next, std::unique_ptr<DecodingState>&) const
{
	//const int minCharCount = 15; // 23 with half-frame number, 15 without
	const int minCharCount = 5;
	const int minQuietZone = 0.2; // Even 1, but it might mess up detection
	//std::cout << "DXFilmEdge::decodePattern()";
	next = FindLeftGuard(next, minCharCount, START_PATTERN_, minQuietZone);
	if (!next.isValid())
		return {};

	int xStart = next.pixelsInFront();

	//std::cout << "Found a valid start pattern." << std::endl;

	// Skip the start pattern (black white black white black
	// The first signal bar is always white
	// (separation between the start pattern and the product number) 
	/*next.shift(5);*/

	//next = next.subView(5);

	next.shift(5);

	//// Detect old type of DX Edge Code (without half-frame number)
	//if (IsRightGuard(next.subView(24), STOP_PATTERN_, minQuietZone)) {
	//	std::cout << "DX Edge code with half-frame number" << std::endl;
	//} else if (IsRightGuard(next.subView(15), STOP_PATTERN_, minQuietZone)) {
	//	std::cout << "DX Edge code without half-frame number" << std::endl;
	//} else {
	//	//std::cout << "Can't guess the DX Edge code type (with or without half frame number)" << std::endl;
	//}
	////next.subView()

	auto signal_begin = next;
	int signal_length = 0;
	while (next.isValid()) {
		signal_length += *next.data() / 8;
		next.shift(1);
	}

	// Check the stop pattern exists and is about to be reached
	next = next.subView(0, 3);
	if (!IsRightGuard(next, STOP_PATTERN_, minQuietZone)) {
		return {};
	}
		


	//std::cout << "Signal length: " << signal_length << std::endl;
	bool signal_with_half_frame_number = true;
	if (signal_length == 23) {
		//std::cout << "DX Edge code with half-frame number" << std::endl;
	} else if (signal_length == 15) {
		//std::cout << "DX Edge code without half-frame number" << std::endl;
		signal_with_half_frame_number = false;
		//return {};
	} else {
		//std::cout << "Unknown type of DX Edge" << std::endl;
		return {};
	}
	next = signal_begin;

	std::vector<bool> signal_data;
	signal_data.reserve(signal_length);
	//auto signal_data_it = signal_data.begin();

	//std::cout << "BEGIN DECODING" << std::endl;

	if (!next.isValid())
		return {};

	 //Blank Separator
	//auto current_bar_width = *next.data() / 8;

	// Signal starts with a white bar (we include the separator just after the start pattern)
	bool current_signal_is_black = false;

	while (next.isValid()) {
		//std::cout << "sum(): " << next.sum() << std::endl;
		//std::cout << "end(): " << next.end() << std::endl;
		//std::cout << "pixelsInFront(): " << next.pixelsInFront() << std::endl;
		//std::cout << "pixelsTillEnd(): " << next.pixelsTillEnd() << std::endl;
		auto current_bar_width = *next.data() / 8;
		while (current_bar_width > 0) {
			if (current_signal_is_black) {
				signal_data.push_back(true);
			} else {
				signal_data.push_back(false);
			}
			--current_bar_width;
		}
		//std::cout << "data: " << current_bar_width << std::endl;
		next.shift(1);
		current_signal_is_black = !current_signal_is_black;
	}
	if (signal_data.at(0) || signal_data.at(8))
			return {};

	if (signal_with_half_frame_number && (signal_data.at(20) || signal_data.at(22)))
		return {};

	if (!signal_with_half_frame_number && (signal_data.at(12) || signal_data.at(14)))
		return {};

	// Check this is not the clock. It can be scanned when the minQuietZone is too small.
	if (minQuietZone <= 1 && signal_with_half_frame_number) {
		std::vector<bool> signal_clock = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0};
		if (signal_clock == signal_data)
			return {};
	}

	// Check parity bit
	auto signal_sum = std::accumulate(signal_data.begin(), signal_data.end()-2, 0);
	//std::cout << "Signal sum: " << signal_sum << std::endl;
	//std::cout << "Parity bit: " << (int)*(signal_data.end() - 2) << std::endl;
	auto parity_bit = *(signal_data.end() - 2);
	if (signal_sum % 2 != parity_bit)
		return {};

	//std::cout << "Signal bitfield: ";
	//for(bool bitfield : signal_data) {
		//std::cout << (int)bitfield;
	//}
	//std::cout << std::endl;
	//std::cout << "Generating Product number..." << std::endl;
	auto product_number = fromVector(signal_data.begin() + 1, signal_data.begin() + 8);
	if (!product_number)
		return {};
	//std::cout << "Generating Generation number..." << std::endl;
	auto generation_number = fromVector(signal_data.begin() + 9, signal_data.begin() + 13);
	//std::cout << "Product number: " << product_number << std::endl;
	//std::cout << "Generation number: " << generation_number << std::endl;
	//if (signal_with_half_frame_number) {
		//auto half_frame_number = fromVector(signal_data.begin() + 13, signal_data.begin() + 20);
		//std::cout << "Frame number: " << half_frame_number / 2;
		//if (half_frame_number % 2)
			//std::cout << "A";
		//std::cout << std::endl;
	//}

	//std::cout << "Suggested DX number: 0" << std::setfill('0') << std::setw(4) << 16 * product_number + generation_number;
	//std::cout << "3" << std::endl;

	//std::cout << "END OF DECODING" << std::endl;
	std::string txt;
	txt.reserve(12);
	txt = std::to_string(product_number) + "-" + std::to_string(generation_number);
	if (signal_with_half_frame_number) {
		auto half_frame_number = fromVector(signal_data.begin() + 13, signal_data.begin() + 20);
		txt += "/" + std::to_string(half_frame_number / 2);
		if (half_frame_number % 2)
			txt += "A";
	}
	//constexpr int weights[] = {1, 2, 4, 7, 0};
	//std::cout << "Total signal sum: " << next.sum() << std::endl;
	//std::cout << "Decoding the signal..." << std::endl;
	//while (next.isValid() && !next.isAtLastBar()) {
	//	txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//	next.shift(1);
	//}
	//std::cout << "Decoded signal: " << txt << std::endl;

	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;

	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;

	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;

	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;

	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;

	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;

	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current index: " << next.index();
	//next.shift(1);
	//std::cout << "current text:" << txt << std::endl;



	//std::cout << "current text:" << txt << std::endl;
	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current text:" << txt << std::endl;
	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current text:" << txt << std::endl;
	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current text:" << txt << std::endl;
	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current text:" << txt << std::endl;
	//txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//std::cout << "current text:" << txt << std::endl;

	//do {
	//	txt += DecodeNarrowWidePattern(next, CHARACTER_ENCODINGS, ALPHABET);
	//	std::cout << "current text:" << txt << std::endl;
	//	if (txt.back() == 0)
	//		return {};
	//} while (!isStartOrStopSymbol(txt.back()));



	/*next = next.subView(4, 10);

	while (next.isValid()) {
		const auto threshold = NarrowWideThreshold(next);
		if (!threshold.isValid())
			break;

		BarAndSpace<int> digits, numWide;
		for (int i = 0; i < 10; ++i) {
			if (next[i] > threshold[i] * 2)
				break;
			numWide[i] += next[i] > threshold[i];
			digits[i] += weights[i / 2] * (next[i] > threshold[i]);
		}

		if (numWide.bar != 2 || numWide.space != 2)
			break;

		for (int i = 0; i < 2; ++i)
			txt.push_back(ToDigit(digits[i] == 11 ? 0 : digits[i]));

		next.skipSymbol();
	}

	next = next.subView(0, 3);*/

	//if (Size(txt) < minCharCount || !next.isValid())
	//	return {};

	//if (Size(txt) < minCharCount)
	//	return {};




	Error error;
	//if (_hints.validateITFCheckSum() && !GTIN::IsCheckDigitValid(txt))
	//	error = ChecksumError();

	// Symbology identifier ISO/IEC 16390:2007 Annex C Table C.1
	// See also GS1 General Specifications 5.1.3 Figure 5.1.3-2
	SymbologyIdentifier symbologyIdentifier = {'I', '0'}; // No check character validation

	// if (_hints.validateITFCheckSum() || (txt.size() == 14 && GTIN::IsCheckDigitValid(txt))) // If no hint test if valid ITF-14
		// symbologyIdentifier.modifier = '1'; // Modulo 10 symbol check character validated and transmitted

	int xStop = next.pixelsTillEnd();
	return Result(txt, rowNumber, xStart, xStop, BarcodeFormat::DXFilmEdge, symbologyIdentifier, error);
}

} // namespace ZXing::OneD
