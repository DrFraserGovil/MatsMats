#pragma once
#include "../biology/MotifMatrix.h"
#include "ScanRecord.h"


using scoreholder = double;

class SequenceScanner
{
	public:
		// std::vector<MotifMatrix> OnTheFly;
		// std::vector<std::vector<MotifMatrix>> Precomputed;
		
		void LoadMotifs(int sequenceCount, int sequenceLength);

		void Scan(Sequence::DNA & dna, ThreadRecords & records);
		size_t size() const;
	private:
		std::vector<std::vector<int>> Precomputers;
		std::vector<int> PrecomputedSizes;
		std::vector<int> Fliers;
		std::vector<MotifMatrix> Motifs;
		int NMotifs;

		//first index groups motifs of the same length (small, < 5)
		//second index is motifs in the group (~100)
		//third index is the dnabits encoding 
		//pair element is the forward and reverse complement of the dnabits
		std::vector<std::vector<std::vector<scoreholder>>> PrecomputedScores;
		void Precompute();
};

bool PrecomputationAllowed(size_t sequenceCount, size_t meanSize, size_t motifLength, int callingID);
