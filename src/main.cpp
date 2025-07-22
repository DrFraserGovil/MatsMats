#include <iostream>

#include "settings/Settings.h"

#include "tools/tools.h"
#include "parallel/parallel.h"
#include "biology/DNASequence.h"
#include "scan/MotifSet.h"
#include "scan/Scanner.h"

int main(int argc, char**argv)
{
	LOG(WARN) << "MATSMATS-ultra version";
	srand(1);
	// srand(time(NULL));
	rand();rand(); //flush the seed

	//interprets the command line arguments, and saves them into the settings object
	Settings.Initialise(argc,argv);

	int Nseq = Settings.System.InputSize;
	int SeqLength = 25;

	SequenceScanner ms;
	ms.LoadMotifs(Nseq,SeqLength);
	ThreadRecords recs(ms.size());

	
	
	std::vector<char> lib = {'a','c','g','t'};
	std::string seq = std::string(SeqLength,' ');
	
	for (int q = 0; q < SeqLength; ++q)
	{
		int r = rand() % 4;
		seq[q] = lib[r];
	}
	// LOG(DEBUG) << seq;
	Sequence::DNA dna(seq);
	int progSize = 16;
	LOG(INFO) << "[" << std::string(progSize,' ') <<"]";
	
	int cTicks = 0;
	double ticksPerSeq = progSize* 1.0/(Nseq-1);
	for (int i = 0; i < Nseq ; ++ i)
	{
		
		for (int q= 0; q < SeqLength; ++q)
		{
			dna.Sequence[q] = rand() % 4;
		}

		int nt = ticksPerSeq * i;
		if (nt > cTicks)
		{
			cTicks = nt;
			LOG(INFO).ErasePrevious();
			LOG(INFO) << "[" << std::string(cTicks,'#') << " " << std::string(progSize - cTicks,' ') <<"]";
		}
		ms.Scan(dna,recs);
	}
	
	



	//Load PWMs from memory

	//Load sequences from memory & pre-analyse are we assuming from file and not piped in? Pipes may be more efficient than saving to disk if preprocessing

	//Then run the analysis

	return 0;
}