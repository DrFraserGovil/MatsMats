#include <iostream>

#include "settings/Settings.h"

#include "tools/tools.h"
#include "parallel/parallel.h"
#include "biology/DNASequence.h"
#include "scan/SequenceScanner.h"

#include "filesystem"
#include "parallel/parallel.h"
#include "scan/fastqReader.h"
#include <mutex>

namespace fs = std::filesystem;

std::string outputName(fs::directory_entry input, const fs::path inputRoot, const fs::path outputRoot)
{
	fs::path relative_path = input.path().lexically_relative(inputRoot);
	fs::path outputFile = outputRoot / relative_path;

	while (outputFile.has_extension())
	{
		outputFile = outputFile.replace_extension("");;
	}
	outputFile = outputFile.replace_extension(".out");

	fs::create_directories(outputFile.parent_path());
	return outputFile.string();
}


void scanSeq()
{

	
	auto pwm = getRecursiveFileList(Settings.Input.PFMDirectory,Settings.Input.PFMRegex);
	auto scanner = SequenceScanner(pwm);

	const fs::path inputRoot(Settings.Input.ReadDirectory.Value());
	const fs::path outputRoot(Settings.Output.OutputDirectory.Value());
	auto fastqFiles = getRecursiveFileList(inputRoot,Settings.Input.ReadRegex);
	// fastqFiles.resize(15);
	
	
	ParallelPool Parallel(Settings.System.ParallelThreads);

	std::atomic<int> globalCount;
	std::mutex lock;

	
	ProgressBar PB(fastqFiles.size(),"Iterating through " +std::to_string(fastqFiles.size()) + " files\n");
	Parallel.For(fastqFiles.size(),[&](int i)
	{
		++globalCount;
		lock.lock();
		PB.Update(globalCount);
		lock.unlock();

		auto extension = fastqFiles[i].path().extension().string();

		auto outname =  outputName(fastqFiles[i],inputRoot,outputRoot);

		std::ofstream outstream(outname);
		if (extension == ".gz")
		{
			gzfastQScan(fastqFiles[i].path().string(),scanner,outstream);
		}
		else
		{
			fastqScan(fastqFiles[i].path().string(),scanner,outstream);
		}
		outstream.close();

	});
	LOG(INFO) << "Scan complete, exiting scope";
}



int main(int argc, char**argv)
{
	//interprets the command line arguments, and saves them into the settings object
	Settings.Initialise(argc,argv);


	// auto pwm = getRecursiveFileList(Settings.Input.PFMDirectory,Settings.Input.PFMRegex);
		
	// int Nseq = Settings.System.InputSize;
	// int SeqLength = 25;
	// auto scanner = SequenceScanner(pwm,Nseq,SeqLength);

	// SuffixTree tree(scanner.Motifs,scanner.Precomputers[2]);

	scanSeq();
	//Load PWMs from memory

	//Load sequences from memory & pre-analyse are we assuming from file and not piped in? Pipes may be more efficient than saving to disk if preprocessing

	//Then run the analysis

	return 0;
}