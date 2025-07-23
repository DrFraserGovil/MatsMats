#include <iostream>

#include "settings/Settings.h"

#include "tools/tools.h"
#include "parallel/parallel.h"
#include "biology/DNASequence.h"
#include "scan/SequenceScanner.h"

#include "filesystem"
int main(int argc, char**argv)
{
	//interprets the command line arguments, and saves them into the settings object
	Settings.Initialise(argc,argv);

	auto fastqFiles = getRecursiveFileList(Settings.Input.ReadDirectory,Settings.Input.ReadRegex);
	auto pwm = getRecursiveFileList(Settings.Input.PFMDirectory,Settings.Input.PFMRegex);
	LOG(INFO) << fastqFiles.size() << " fastq files";
	LOG(INFO) << pwm.size() << " pwm files";

	
	int Nseq = Settings.System.InputSize;
	int SeqLength = 25;
	auto scanner = SequenceScanner(pwm,Nseq,SeqLength);

	// SequenceScanner ms;
	// ms.LoadMotifs(Nseq,SeqLength);

	// srand(1);
	
	std::vector<char> lib = {'a','c','g','t'};
	std::string seq = std::string(SeqLength,' ');
	
	for (int q = 0; q < SeqLength; ++q)
	{
		int r = rand() % 4;
		seq[q] = lib[r];
	}
	// LOG(DEBUG) << seq;
	Sequence::DNA dna(seq);
	
	Record best;
	Record bestbest;
	Sequence::DNA bestDNA(seq);
	double bestScore = -99999;
	Progress::For(0,Nseq,[&](int i)
	{
		for (int q= 0; q < SeqLength; ++q)
		{
			int r=  rand() % 4;
			dna.Sequence[q] = r;
			dna.SequenceString[q] = lib[r]; 
			// seq[q] = lib[rand() % 4];
		}
		// dna.NewSequence(seq);
		scanner.Scan(dna,best);
		if (best.Score > bestScore)
		{
			bestbest= best;
			bestScore = best.Score;
			bestDNA = dna;
		}
		
	},"Scanning Sequences\n");


	



	//Load PWMs from memory

	//Load sequences from memory & pre-analyse are we assuming from file and not piped in? Pipes may be more efficient than saving to disk if preprocessing

	//Then run the analysis

	return 0;
}