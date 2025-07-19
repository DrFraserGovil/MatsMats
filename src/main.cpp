#include <iostream>

#include "tools/tools.h"
#include "parallel/parallel.h"


void incrementFunction(int i,std::atomic<int>* counter)
{
	(*counter)++;
	LOG(INFO) << *counter;
}
int main(int argc, char**argv)
{
	for (int nWorkers = 1; nWorkers < 10; ++nWorkers)
	{
		auto Parallel = ParallelPool(nWorkers);
		std::atomic<int> counter =0;
		int nIncrements = 100;
		Parallel.For(nIncrements,incrementFunction,&counter);
		LOG(ERROR) << nWorkers << " " << counter << " " << nIncrements;
	}

	return 0;
}