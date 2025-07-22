#include "MotifMatrix.h"
#include <iomanip>
#include <cmath>

MotifMatrix::MotifMatrix(std::string filepath,int id) : ID(id)
{
	//placeholder random initialisation
	int L = 10;//4 + rand() % 10;
	LogOdds = std::vector<std::vector<double>>(L,std::vector<double>(4,0.0));

	for (int i = 0; i < L; ++i)
	{
		double n = 0;
		for (int j = 0; j < 4; ++j)
		{
			LogOdds[i][j] = -4*(rand()*1.0/RAND_MAX -1);
			n += std::exp(LogOdds[i][j]);
		}
		n = std::log(n);
		for (int j = 0; j < 4; ++j)
		{
			LogOdds[i][j] -= n + std::log(0.25);
		}
	}
	MotifLength = LogOdds.size();
}
size_t MotifMatrix::size() const
{
	return MotifLength;
}




const std::vector<int> rcMap = {3,2,1,0}; //indexed so that rcMap[base] = rc_base
std::pair<double,double> MotifMatrix::Score(Sequence::DNA & sequence, int idx) const
{
	double forwardScore = 0;
	double rcScore = 0;
	const int L = MotifLength;
	for (int i = 0; i < L; ++i)
	{
		int base = sequence.Sequence[idx+i];
		forwardScore += LogOdds[i][base];
		rcScore += LogOdds[L-i-1][rcMap[base]];
	}
	return {forwardScore/L,rcScore/L};
}

std::string MotifMatrix::ToString() const
{
	std::stringstream s;
	for (int i =0; i < 4; ++i)
	{
		for (int j = 0; j < MotifLength; ++j)
		{
			s << LogOdds[j][i];
			if (j < MotifLength -1)
			{
				s << " ";
			}
		}
		if (i < 3)
		{
			s << "\n";
		}
	}
	return s.str();
}

const std::vector<std::vector<double>> & MotifMatrix::ReferenceScores() const
{
	return LogOdds;
}