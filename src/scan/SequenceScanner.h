#pragma once
#include "../biology/MotifMatrix.h"
#include "ScanRecord.h"
#include <filesystem>

using fs_path = std::filesystem::directory_entry;

class SequenceScanner
{
	public:
		// std::vector<MotifMatrix> OnTheFly;
		// std::vector<std::vector<MotifMatrix>> Precomputed;
		SequenceScanner(std::vector<fs_path> motifPaths, int sequenceCount=Settings.Input.EstimatedReadCount,int sequenceLength=Settings.Input.EstimatedReadLength);
		
		void Scan(Sequence::DNA & dna, Record & record);
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
		void InitialiseMotifs(int sequenceCount, int sequenceLength);
		void Precompute();
};

bool PrecomputationAllowed(size_t sequenceCount, size_t meanSize, size_t motifLength, int callingID);
