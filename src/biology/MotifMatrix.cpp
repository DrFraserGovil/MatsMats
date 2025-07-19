#include "MotifMatrix.h"
#include "iomanip"
MotifMatrix::MotifMatrix(std::string filepath,int id)
{
	ID = id;
	//placeholder random initialisation
	int L = 12;
	LogOdds = std::vector<std::vector<double>>(L,std::vector<double>(4,0.0));

	for (int i = 0; i < L; ++i)
	{
		double n = 0;
		for (int j = 0; j < 4; ++j)
		{
			LogOdds[i][j] = -4*(rand()*1.0/RAND_MAX -1);
			n += exp(LogOdds[i][j]);
		}
		n = log(n);
		for (int j = 0; j < 4; ++j)
		{
			LogOdds[i][j] -= n + log(0.25);
		}
	}
}

void MotifMatrix::Initialise(size_t sequenceCount, size_t meanSize)
{
	int MotifLength = LogOdds.size();

	//Check if precomputing would actually speed things up
	long long int PrecomputeSize = pow(4,MotifLength); // number of k-mers required for precomputing
	int ExpectedKmerCount = sequenceCount *  std::max(1,(int)(meanSize - LogOdds.size())); 
	bool isFaster = (ExpectedKmerCount * 1.0/PrecomputeSize > 1);

	//Check if precomputing possible within encoding limitations
	bool smallEnoughForEncoding = (Sequence::MaximumEncodingLength() >= MotifLength);

	//Check if memory footprint would be ludicrous
	//footprint factor estimates global memory from this single factor, chosen by user (default = 100)
	double expectedGlobalFootprint = (PrecomputeSize * sizeof(double)) *1.0/ pow(1024,3) * Settings.System.FootprintFactor;
	bool fitsInMemory = (expectedGlobalFootprint <= Settings.System.MemoryLimit);
	
	PrecomputeMode = (isFaster) && (smallEnoughForEncoding) && (fitsInMemory);


	if (Settings.System.Verbosity >= LogLevel::DEBUG)
	{
		//output a nicely formatted table reasoning why precomputation was initialised (or not)
		std::stringstream buffer;
		buffer << "Initialising Motif " << ID;
		buffer << "\n   Length:   " << std::setw(6) << MotifLength;
		buffer << "\tPC-Kmers: " << std::setw(12)<< PrecomputeSize;
		buffer << "\tEst-Kmers: " << std::setw(10) << ExpectedKmerCount;
		buffer << "\tFootprint: " << std::setw(10) << expectedGlobalFootprint << "GiB";
		buffer << "\n   Speedgain:" << std::setw(6) << MakeString(isFaster);
		buffer << "\tEncoder valid:" << std::setw(8) << MakeString(smallEnoughForEncoding);
		buffer << "\tIn memory: " << std::setw(10) << MakeString(fitsInMemory);
		buffer << "\n   ";
		if (PrecomputeMode)
		{
			buffer << "PRECOMPUTE";
		}
		else
		{
			buffer << "ON-THE-FLY";
		}
		LOG(DEBUG) << buffer.str();
	}

	if (PrecomputeMode)
	{
		PrecomputeScores();
	}
}

void MotifMatrix::PrecomputeScores()
{
	int L = LogOdds.size();
	dnabits nCodes = pow(4,L);
	PrecomputedScores.resize(nCodes);
	double bs = -99999;
	std::string bestSeq;
	for (dnabits code = 0; code < nCodes; ++code)
	{
		double score = 0;
		dnabits decoder = code;
		for (int i = 0; i < L; ++i)
		{
			int base = decoder & 3;
			decoder = decoder >> 2;
			score += LogOdds[L-i-1][base];
		}

		PrecomputedScores[code] = score;
		if (score > bs)
		{
			bs = score;
			bestSeq = Sequence::Decode(code,L);
		}
	}
	LOG(DEBUG) << "Precomputation completed. Best matching motif is " << bestSeq;
}

double MotifMatrix::BestScore(Sequence::DNA & input)
{
	// if (PrecomputeMode)
	// {
	// 	input.ResetBitfield();
	// }
	// return 0;
	//iterate over all subsequences 
}

// double MotifMatrix::Score()
// {
// 	double score = 0;

// }