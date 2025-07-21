#pragma once
#include <string>
#include <vector>
#include "../tools/Log.h"
enum Direction
{
	Forward,
	Backward
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

struct Record
{
	double Score;
	int Position;
	Direction Direction;
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
	// private:

};