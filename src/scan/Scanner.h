#pragma once

#include "../biology/DNASequence.h"
#include "MotifSet.h"

#include "ScanRecord.h"



void checkScore(Record & runningBest, double score, int pos, Direction direction,bool forceWin)
{
	if (score > runningBest.Score || forceWin)
	{
		runningBest.Score = score;
		runningBest.Position = pos;
		runningBest.Direction = direction;
		runningBest.Hits = 1;
	}
	else if (abs((score - runningBest.Score)) < 1e-8)
	{
		runningBest.Hits +=1;
	}
}



void FlyingScan(Sequence::DNA & dna, const std::vector<MotifMatrix> & motifs, std::vector<Record> & records)
{
	for (int j = 0; j < motifs.size(); ++j)
	{
		
	}
}

// void Scan(Sequence::DNA & dna, const MotifSet & motifs)
// {
// 	std::vector<Record> motifRecords(motifs.size());
// 	double globalBestID = -1;
// 	double globalBestScore = 0;
// 	/*
// 		There's probably some optimisation to be had if we scanned each position with each motif (rather than do a full scan per motif)
// 		This would avoid recalculating the bitfields. However, I imaging that this is a subdominant factor, so go for the simple approach first
// 	*/
// 	// for (int i = 0; i < motifs.size(); ++i)
// 	// {
// 	// 	// LOG(INFO) << i;
// 	// 	Record runningBest;

// 	// 	auto & motif = motifs[i];
// 	// 	int scanSize= dna.Length - motif.size() +1; 
// 	// 	dna.SetBitfield(motif.size());
// 	// 	for (int j =0; j < scanSize; ++j)
// 	// 	{
// 	// 		// motif.Score(dna,j);
// 	// 		// LOG(INFO) << j << " " << Sequence::Decode(dna.GetBitfield(),motif.size()) << " " <<dna.Score << " " << Sequence::Decode(dna.GetRCBitfield(),motif.size()) << dna.RCScore;
			
// 	// 		// checkScore(runningBest,dna.Score,j,Direction::Forward,j==0);
// 	// 		// checkScore(runningBest,dna.RCScore,j,Direction::Backward,false);

			
// 	// 		// if (j < scanSize - 1) // on all but the last step
// 	// 		// {
// 	// 		// 	dna.StepBitfield();
// 	// 		// }
// 	// 	}

// 	// 	// motifRecords[i] = runningBest;
// 	// 	if (runningBest.Score > globalBestScore || i == 0)
// 	// 	{
// 	// 		globalBestScore = runningBest.Score;
// 	// 		globalBestID = i;
// 	// 	}
// 	// }

// 	// auto best = motifRecords[globalBestID];
// 	// LOG(INFO) << "Motif " << globalBestID << " was the winner, with a score " << globalBestScore << " " << best.Direction << " " << best.Hits;
	
// 	// std::string place = "[string]";
	
// 	// dna.FileString = std::format("{} {:3} {:3} {} {:2} {:7.6f}",place,globalBestID,best.Position,directionString(best.Direction),best.Hits,best.Score);
// 	// if (Settings.Output.FullSummary)
// 	// {
// 	// 	for (int i = 0; i < motifRecords.size(); ++i)
// 	// 	{
// 	// 		auto & record = motifRecords[i];
// 	// 		std::format_to(std::back_inserter(dna.FileString)," | {:3} {} {:.6f}",record.Position,directionString(record.Direction),record.Score);
// 	// 	}
// 	// }
// 	// if (id == 99 )
// 	// {
// 	// 	LOG(INFO) << id << " " << dna.FileString;
// 	// }
// 	// ++id;
// };