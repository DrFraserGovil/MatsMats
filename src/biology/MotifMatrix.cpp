#include "MotifMatrix.h"
#include <iomanip>
#include <cmath>
#include "../tools/fileparser.h"
MotifMatrix::MotifMatrix(std::string filepath,int id) : ID(id)
{
	LOG(INFO) << "Loading " << filepath;
	int L = -1;
	int lineno = 1;
	int expectedLineNo = 4;
	forSplitLineIn(filepath,"\t",[&](auto line)
	{
		if (lineno > expectedLineNo)
		{
			LOG(ERROR) << "This PFM file has more than " << expectedLineNo << " lines, and is therefore not a proper DNA-PFM";
			throw std::runtime_error("Improperly formated PFM");
		}
		if (L < 0)
		{
			L = line.size();
			LogOdds = std::vector<std::vector<double>>(L,std::vector<double>(4,0.0));		
		}
		else
		{
			if (line.size() != L)
			{
				LOG(ERROR) << "PFM " << filepath << "has an irregular length";
				throw std::runtime_error("Improperly formated PFM");
			}
		}
		for (int i = 0; i < L; ++i)
		{
			LogOdds[i][lineno-1] = convert<double>(line[i]);
		}
		++lineno;
		// LOG(INFO) << MakeString(line);

	});
	if (lineno -1 != expectedLineNo)
	{
		LOG(ERROR) << "This PFM file had only " << lineno -1 << " lines, and is therefore not a proper DNA-PFM";
		throw std::runtime_error("Improperly formated PFM");
	}
	//placeholder random initialisation
	// int L = 10 + rand() % 8;
	// LogOdds = std::vector<std::vector<double>>(L,std::vector<double>(4,0.0));

	for (int i = 0; i < L; ++i)
	{
		double n = 0;
		for (int j = 0; j < 4; ++j)
		{
			LogOdds[i][j] = std::log(LogOdds[i][j]);
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