#pragma once
#include <string_view>
#include <vector>
#include "../tools/Log.h"
typedef uint_fast32_t dnabits;


namespace Sequence
{
	const int LogAlphabetSize = 2;
	const int BitHackExtractor = 3;
	class DNA
	{
		public:
			DNA(std::string_view sequencestring);

			int Length;
			std::vector<unsigned char> Sequence; // array of 0/1/2/3 for A/C/T/G. Unsigned char is 1byte long, so nice and efficient

			//bitfield manipulators
			void SetBitfield(size_t motifSize, int startidx=0);
			dnabits GetBitfield();
			dnabits GetRCBitfield();
			void StepBitfield();
			
			void NewSequence(std::string_view sequence);
			static int MaximumEncodedLength();

			// int SubstringHead;
			// double Score;
			// double RCScore;
			std::string FileString;
			std::string SequenceString;
		private:

			
			dnabits Bitfield;
			size_t FieldRightIndex;
			int CurrentMotifSize;
			
			dnabits Mask;

	};


	/*! @brief Calculates the maximum length of a subsequence that can be encoded in a dnabits field
		@details dnabits is defined to be uint_fast32, which means that it is always at least 32 bits (but may be more, system dependent). This means the maximum encoding length is not always certain. If more encoding space is required, upgrade dnabits to uint_fast64 to double the capacity.
		@returns A length (in bases) above which dnabits representations are no longer valid
	*/

	constexpr size_t MaximumEncodingLength()
	{
		constexpr int DNA_BASES_PER_BYTE = 4; // 8 bits/byte / 2 bits/base = 4 bases/byte
		return sizeof(dnabits) * DNA_BASES_PER_BYTE;
	}


	std::string Decode(dnabits code, size_t length);
}