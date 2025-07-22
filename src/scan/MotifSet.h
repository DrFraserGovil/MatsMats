#pragma once
#include "../biology/MotifMatrix.h"
#include "ScanRecord.h"


using scoreholder = short int;


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
		//second index is the dnabits encoding 
		std::vector<std::vector<PrecomputeElement>> PrecomputedScores;
		void Precompute();
};

bool PrecomputationAllowed(size_t sequenceCount, size_t meanSize, size_t motifLength, int callingID);
