#include <iostream>

#include "settings/Settings.h"

#include "tools/tools.h"
#include "parallel/parallel.h"
#include "biology/DNASequence.h"
#include "biology/MotifMatrix.h"

int main(int argc, char**argv)
{
	//interprets the command line arguments, and saves them into the settings object
	Settings.Initialise(argc,argv);

	std::vector<MotifMatrix> ms;
	for (int id = 0; id < 100; ++id)
	{
		auto motif = MotifMatrix("test",id);
		motif.Initialise(1e7,50);
		ms.push_back(motif);
	}
	LOG(INFO) << ms.size();

	//Load PWMs from memory

	//Load sequences from memory & pre-analyse are we assuming from file and not piped in? Pipes may be more efficient than saving to disk if preprocessing

	//Then run the analysis

	return 0;
}