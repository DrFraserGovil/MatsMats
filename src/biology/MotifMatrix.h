#pragma once
#include <vector>
#include "../tools/tools.h"
#include "DNASequence.h"
class MotifMatrix
{
	public:
		int Length;

		double BestScore(DNASequence & input);
		bool PrecomputeMode;
	private:
		std::vector<std::vector<double>> Matrix;

		std::vector<double> ScoreList;
};