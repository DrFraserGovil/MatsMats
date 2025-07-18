#include "DNASequence.h"


DNASequence::DNASequence(std::string_view sequencestring): Length(sequencestring.size())
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
}

dnabits DNASequence::ResetBitfield(size_t motifSize)
{
	 if (motifSize > Length)
	 {
		LOG(ERROR) << "Motif size (" << motifSize << ") exceeds sequence length (" << Length <<"). MATSMATS does not allow this.";
        throw std::out_of_range("Motif size exceeds sequence length.");
    }

	CurrentMotifSize = motifSize;
	FieldRightIndex = motifSize-1;
	Bitfield = 0;
	for (int i = 0; i < motifSize; ++i)
	{
		Bitfield= (Bitfield << LogAlphabetSize) + Sequence[i];
	}
	Mask = (static_cast<dnabits>(1) << (LogAlphabetSize * motifSize)) - 1; //mask must be of correct type, take no prisoners!
	return Bitfield;
}
dnabits DNASequence::GetBitfield()
{
	return Bitfield;
}

dnabits DNASequence::NextBitfield()
{
	//performs no checks of its own (for speed), so caller responsibility to make sure FieldRightIndex is never out of bounds 
	++FieldRightIndex;
	Bitfield = ((Bitfield << LogAlphabetSize) + Sequence[FieldRightIndex]) & Mask;
	return Bitfield;
}