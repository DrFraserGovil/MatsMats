#pragma once
#include <string>
#include <vector>
#include "../tools/Log.h"
enum Direction
{
	Forward,
	Backward,
	Uninitialised,
};

std::string inline directionString(Direction dir)
{
	switch(dir)
	{
		case Direction::Forward:
			return  "+";
			break;
		case Direction::Backward:
			return  "-";
			break;
		default:
			return "?";
			break;
	}
}

struct PrecomputeElement
{
	double Score;
	Direction Strand;
	int MotifID;

	PrecomputeElement(){Strand = Uninitialised;};

	void CheckElement(double fscore, double rcscore,int motif)
	{
		if (fscore > Score || Strand == Uninitialised)
		{
			Score = fscore;
			Strand = Forward;
			MotifID = motif;
		}
		if (rcscore > Score)
		{
			Score = rcscore;
			Strand = Backward;
			MotifID = motif;
		}
	};
};

struct Record
{
	double Score;
	int Position;
	Direction Strand;
	int Hits;
};

class ThreadRecords
{
	public:
		ThreadRecords(int nMotifs);
		int BestID;
		double BestScore;
		std::vector<Record> Records;

		void CheckRecords(int motifID, double score, int pos, Direction dir, bool forceWin);
		void NewScan();

};