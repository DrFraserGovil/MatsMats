#pragma once
#include <vector>
#include "../tools/tools.h"
#include "DNASequence.h"
#include "../settings/Settings.h"
class MotifMatrix
{
	public:
		MotifMatrix() : ID(-1){};
		MotifMatrix(std::string path,int id);
		// double BestScore(Sequence::DNA & input);
		// void Initialise(size_t SequenceCount, size_t MeanLength);
		size_t size() const;
		std::string ToString() const;
		
		std::vector<std::vector<double>> LogOdds;
		const int ID;
		std::pair<double,double> Score(Sequence::DNA & sequence, int idx) const;
		const std::vector<std::vector<double>> & ReferenceScores() const;
	private:
		size_t MotifLength;
		// void PrecomputeScores();
		// std::vector<double> PrecomputedScores;
};



/*
	PLAN -- simplify this into a basic class (for on the fly) and then a holder class for an aggregate precomputed scores with motif lengths sizes. The on-the-fly class is lightweight so can be duplicated.


*/