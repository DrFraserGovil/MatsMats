#include "ScanRecord.h"



void Record::CheckRecords(int motifID, double score, int pos, Direction dir, bool forceWin)
{
	if (score > Score || forceWin)
	{
		MotifID = motifID;
		Score = score;
		Position = pos;
		Strand = dir;
		Hits = 1;
	}
	else if (abs(score - Score) < 1e-8)
	{
		++Hits;
	}

}

