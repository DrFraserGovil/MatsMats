#pragma once
#include <vector>
#include "../tools/tools.h"
#include "DNASequence.h"
#include "../settings/Settings.h"
class MotifMatrix
{
	public:
		MotifMatrix(std::string path,int id);
		double BestScore(Sequence::DNA & input);
		void Initialise(size_t SequenceCount, size_t MeanLength);
		size_t size();
	private:
		bool PrecomputeMode;
		void PrecomputeScores();
		std::vector<std::vector<double>> LogOdds;
		std::vector<double> PrecomputedScores;
		int ID;
		double Score(dnabits subsequence);
};