#include "MotifSet.h"
#include <iomanip>
#include <cmath>
#include "../tools/formatter.h"
// #include <format>

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


void SequenceScanner::LoadMotifs(int sequenceCount, int sequenceLength)
{
	Motifs.resize(0);
	for (int i = 0; i < 200; ++i)
	{
		//construct directly into the holder array for efficiency
		Motifs.emplace_back("test",i);

		auto & newMotif = Motifs[i];

		//now determine if the motif will be precomputed or on-the-fly
		if (PrecomputationAllowed(sequenceCount,sequenceLength,newMotif.size(),newMotif.ID))
		{
			//if precomputed, group it with all the other precomputation grids with the same motif length

			//this is just a manually written 'if L in vec, insert at position, otherwise, create new entry'.
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

	//finalise the initialisation - logging and the precomputation
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
	int completed = 0;
	PrecomputedScores.resize(Precomputers.size());
	
	ProgressBar PB(nPrecompute,"Precomputing Score Tables\n");
	for (int i = 0; i < PrecomputedSizes.size(); ++i)
	{
		//all motifs in this section have the same motif length
		int L = Motifs[Precomputers[i][0]].size();
	
		//there are 4^L L-mers, and we're going to iterate over all of them
		//This is why it's important to check that this is feasible! (See: PrecomputationAllowed())
		dnabits nCodes = std::pow(4,L);	
		PrecomputedScores[i].resize(nCodes);

		for (int j = 0; j < Precomputers[i].size(); ++j)
		{
			PB.Update(completed);	
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

				//this is where the magic happens. We compute both the forward and rc score, and then check them against the scores achieved by *all of the motifs in the set*.
				//We then store the winner. We then only need to do a single lookup for each subsequence to learn the best-scoring motif, and the best-scoring direction
				PrecomputedScores[i][code].CheckElement(fscore/L,rcscore/L,Precomputers[i][j]);	
			}
			++completed;
		}
	}
}


void SequenceScanner::Scan(Sequence::DNA & dna, Record & best)
{
	//it's important that firstCheck is only active once, because it force-resets the record.
	//Need to be careful on scans with mixed fly/precomp sets!

	bool firstCheck = true;
	//do the on the fly motifs first
	for (int i = 0; i < Fliers.size(); ++i)
	{
		auto & motif = Motifs[Fliers[i]];
		// auto & record = records[Fliers[i]];
		int scanSize= dna.Length - motif.size() +1; 
		for (int j = 0; j < scanSize; ++j)
		{
			
			auto[fscore,rcscore] = motif.Score(dna,j);
			best.CheckRecords(Fliers[i],fscore,j,Direction::Forward,firstCheck);
			firstCheck = false;
			best.CheckRecords(Fliers[i],rcscore,j,Direction::Backward,false);
		}
	}


	//then do the precomputes
	for (int i = 0; i < Precomputers.size(); ++i)
	{
		int L = PrecomputedSizes[i];
		dna.SetBitfield(L,0);
		int scanSize= dna.Length - L +1; 
		for (int j = 0; j < scanSize; ++j)
		{
			dnabits pos = dna.GetBitfield();

			auto pre = PrecomputedScores[i][pos];
			best.CheckRecords(pre.MotifID,pre.Score,j,pre.Strand,firstCheck);
			firstCheck = false;

			if (j < scanSize - 1) // on all but the last step
			{
				dna.StepBitfield();
			}
		}
	}
	int L = Motifs[best.MotifID].size();
	format_to_string(dna.FileString,"%s %d %d %d %s %d %f",dna.SequenceString,best.MotifID,best.Position,best.Position+L, directionString(best.Strand),best.Hits,best.Score);

}

