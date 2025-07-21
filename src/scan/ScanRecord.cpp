#include "ScanRecord.h"

ThreadRecords::ThreadRecords(int nMotifs)
{
	Records.resize(nMotifs);
}


void ThreadRecords::CheckRecords(int motifID, double score, int pos, Direction dir, bool forceWin)
{
	auto & runningBest = Records[motifID];
	if (score > runningBest.Score || forceWin)
	{
		runningBest.Score = score;
		runningBest.Position = pos;
		runningBest.Strand = dir;
		runningBest.Hits = 1;
	}
	else if (abs((score - runningBest.Score)) < 1e-8)
	{
		runningBest.Hits +=1;
	}

	if (BestID == -1 || score > BestScore)
	{
		BestID = motifID;
		BestScore = score;
	}
}


void ThreadRecords::NewScan()
{
	BestID = -1;
	BestScore = 0;
}