#include "DNASequence.h"

namespace Sequence
{
	//used for the bitfields to efficiently encode dna strings
	

	DNA::DNA(std::string_view sequencestring): Length(sequencestring.size())
	{
		Sequence.reserve(sequencestring.length()); // Pre-allocate memory for efficiency
		for (char c : sequencestring)
		{
			switch (c) 
			{
				case 'A': case 'a': Sequence.push_back(0); break;
				case 'C': case 'c': Sequence.push_back(1); break;
				case 'G': case 'g': Sequence.push_back(2); break;
				case 'T': case 't': Sequence.push_back(3); break;
				default:
					LOG(ERROR) << "Character " << c << " in sequence " << sequencestring << " is not a part of the DNA library";
					throw std::runtime_error("Invalid DNA sequence");
					break;
			}
		}
		// Score = 0;
		// RCScore = 0;
	}

	void DNA::SetBitfield(size_t motifSize, int startIdx )
	{
		if (motifSize > Length)
		{
			LOG(ERROR) << "Motif size (" << motifSize << ") exceeds sequence length (" << Length <<"). MATSMATS does not allow this.";
			throw std::out_of_range("Motif size exceeds sequence length.");
		}

		CurrentMotifSize = motifSize;
		FieldRightIndex = startIdx + motifSize - 1;
		Bitfield = 0;
		for (int i = 0; i < motifSize; ++i)
		{
			Bitfield= (Bitfield << LogAlphabetSize) + Sequence[i+startIdx];
		}
		Mask = (static_cast<dnabits>(1) << (LogAlphabetSize * motifSize)) - 1; //mask must be of correct type, take no prisoners!
	}
	dnabits DNA::GetBitfield()
	{
		return Bitfield;
	}

	dnabits DNA::GetRCBitfield()
	{
		// return 0;//
		dnabits rc = 0;
		dnabits orig = Bitfield;
		
		for (int i = 0; i < CurrentMotifSize; ++i)
		{
			rc = (rc <<LogAlphabetSize) + (orig & BitHackExtractor) ^ BitHackExtractor; //bit hacking!
			orig = orig >> LogAlphabetSize;
		}
		return rc;
	}

	void DNA::StepBitfield()
	{
		//performs no checks of its own (for speed), so caller responsibility to make sure FieldRightIndex is never out of bounds 
		++FieldRightIndex;
		Bitfield = ((Bitfield << LogAlphabetSize) + Sequence[FieldRightIndex]) & Mask;
	}


	std::vector<char> decodermap = {'A','C','G','T'};
	std::string Decode(dnabits code,size_t length)
	{
		std::string out(length,' ');
		for (int i = 0; i < length; ++i)
		{
			char base = decodermap[code &BitHackExtractor];
			code = code >> LogAlphabetSize;
			out[length-i-1] = base;
		}
		return out;
	}
}