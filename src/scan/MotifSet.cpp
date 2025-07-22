#include "MotifSet.h"
#include <iomanip>
#include <cmath>
// #include <format>

void SequenceScanner::LoadMotifs(int sequenceCount, int sequenceLength)
{
	Motifs.resize(0);
	for (int i = 0; i < 100; ++i)
	{
		//construct directly into the holder array for efficiency
		Motifs.emplace_back("test",i);

		auto & newMotif = Motifs[i];
		if (PrecomputationAllowed(sequenceCount,sequenceLength,newMotif.size(),newMotif.ID))
		{
			int L = Motifs[i].size();
			bool found = false;
			for (int j = 0; j < Precomputers.size(); ++j)
			{
				if (PrecomputedSizes[j] == L)
				{
					Precomputers[j].push_back(i);
					found = true;
					break;
					
				}
			}
			if (!found)
			{
				Precomputers.push_back({i});
				PrecomputedSizes.push_back(L);
			}
		}
		else
		{
			Fliers.push_back(i);
		}

	}

	
	NMotifs = Motifs.size();

	if (Fliers.size() > 0)
	{
		LOG(INFO) << Fliers.size() << " on-the-fly arrays found";
	}
	if (PrecomputedSizes.size() > 0)
	{
		LOG(INFO) << NMotifs -  Fliers.size() << " precomputed arrays found:";
		for (int j = 0; j < Precomputers.size(); ++j)
		{
			LOG(DEBUG) << "    " << Precomputers[j].size() << " motifs of size " << PrecomputedSizes[j]; 
		}
		Precompute();
	}
}

size_t SequenceScanner::size() const
{
	return NMotifs;
}

void SequenceScanner::Precompute()
{
	int nPrecompute = NMotifs - Fliers.size();
	int completed = 1;
	double marksPerSeq = 16.0/nPrecompute;
	PrecomputedScores.resize(Precomputers.size());
	
	ProgressBar PB(nPrecompute,"Precomputing Score Tables\n");
	for (int i = 0; i < PrecomputedSizes.size(); ++i)
	{
		
		//all motifs in this section have the same motif length
		int L = Motifs[Precomputers[i][0]].size();
		dnabits nCodes = std::pow(4,L);
		
		PrecomputedScores[i].resize(nCodes);

		for (int j = 0; j < Precomputers[i].size(); ++j)
		{

			auto scores = Motifs[Precomputers[i][j]].ReferenceScores();//const reference to the internal log-odds of the relevant motif
			
			for (dnabits code = 0; code < nCodes; ++code)
			{
				double fscore = 0;
				double rcscore = 0;
				dnabits decoder = code;
				for (int i = 0; i < L; ++i)
				{
					int base = decoder & Sequence::BitHackExtractor;
					int rcbase = base ^ Sequence::BitHackExtractor;
					decoder = decoder >> Sequence::LogAlphabetSize;
					fscore += scores[L-i-1][base];
					rcscore += scores[i][rcbase];
				}

				PrecomputedScores[i][code].CheckElement(fscore,rcscore,Precomputers[i][j]);
				
				
			}
			PB.Update(completed);
			++completed;
			// if (completed <= nPrecompute)
			// {
			// 	LOG(INFO).ErasePrevious();
			// }
		}
	}
}

void SequenceScanner::Scan(Sequence::DNA & dna, ThreadRecords & records)
{
	//do the on the fly motifs first
	records.NewScan();
	for (int i = 0; i < Fliers.size(); ++i)
	{
		auto & motif = Motifs[Fliers[i]];
		// auto & record = records[Fliers[i]];
		int scanSize= dna.Length - motif.size() +1; 
		for (int j = 0; j < scanSize; ++j)
		{
			
			auto[fscore,rcscore] = motif.Score(dna,j);
			records.CheckRecords(Fliers[i],fscore,j,Direction::Forward,j==0);
			records.CheckRecords(Fliers[i],rcscore,j,Direction::Backward,false);
		}
	}


	//then do the precomputes
	for (int i = 0; i < Precomputers.size(); ++i)
	{
		auto & local = PrecomputedScores[i];
		int L = PrecomputedSizes[i];
		dna.SetBitfield(L,0);
		int scanSize= dna.Length - L +1; 
		int nMotifs = Precomputers[i].size();
		for (int j = 0; j < scanSize; ++j)
		{
			dnabits pos = dna.GetBitfield();

			auto best = PrecomputedScores[i][pos];
			records.CheckRecords(best.MotifID,best.Score,j,best.Strand,j==0);

			if (j < scanSize - 1) // on all but the last step
			{
				dna.StepBitfield();
			}
		}
	}
}

bool PrecomputationAllowed(size_t sequenceCount, size_t meanSize, size_t motifLength, int callingID)
{

	//Check if precomputing would actually speed things up
	long long int PrecomputeSize = pow(4,motifLength); // number of k-mers required for precomputing
	long long int ExpectedKmerCount = sequenceCount *  std::max(1,(int)(meanSize - motifLength)); 
	bool isFaster = (ExpectedKmerCount * 1.0/PrecomputeSize > 1);

	//Check if precomputing possible within encoding limitations
	bool smallEnoughForEncoding = (Sequence::MaximumEncodingLength() >= motifLength);

	//Check if memory footprint would be ludicrous
	//footprint factor estimates global memory from this single factor, chosen by user (default = 100)
	double expectedGlobalFootprint = (PrecomputeSize * sizeof(double)) *1.0/ pow(1024,3) * Settings.System.FootprintFactor;
	bool fitsInMemory = (expectedGlobalFootprint <= Settings.System.MemoryLimit);
	
	bool verdict = (isFaster) && (smallEnoughForEncoding) && (fitsInMemory) && !(Settings.System.DisablePrecompute.Value());

	if (Settings.System.Verbosity >= LogLevel::DEBUG)
	{
		//output a nicely formatted table reasoning why precomputation was initialised (or not)
		std::stringstream buffer;
		buffer << "Analysing Motif " << callingID;
		buffer << "\n   Length:   " << std::setw(6) << motifLength;
		buffer << "\tPC-Kmers: " << std::setw(12)<< PrecomputeSize;
		buffer << "\tEst-Kmers: " << std::setw(10) << ExpectedKmerCount;
		buffer << "\tFootprint: " << std::setw(10) << expectedGlobalFootprint << "GiB";
		buffer << "\n   Speedgain:" << std::setw(6) << MakeString(isFaster);
		buffer << "\tEncoder valid:" << std::setw(8) << MakeString(smallEnoughForEncoding);
		buffer << "\tIn memory: " << std::setw(10) << MakeString(fitsInMemory);
		buffer << "\tUser Disabled:" << std::setw(10) << MakeString(Settings.System.DisablePrecompute.Value()); 
		buffer << "\n   Verdict:   ";
		if (verdict)
		{
			buffer << "PRECOMPUTE";
		}
		else
		{
			buffer << "ON-THE-FLY";
		}
		LOG(DEBUG) << buffer.str();
	}
	return verdict;
}