#pragma once
#include <string_view>
#include <vector>
#include "../tools/Log.h"
typedef uint_fast32_t dnabits;

class DNASequence
{
	public:
		DNASequence(std::string_view sequencestring);

		const int Length;
		dnabits ResetBitfield(size_t motifSize);
		dnabits GetBitfield();
		dnabits NextBitfield();
		std::vector<unsigned char> Sequence; // array of 0/1/2/3 for A/C/T/G. Unsigned char is 1byte long, so nice and efficient
	private:

		dnabits Bitfield;
		size_t FieldRightIndex;
		int CurrentMotifSize;
		const int LogAlphabetSize = 2;
		dnabits Mask;
};