#include <iostream>

#include "settings/Settings.h"

#include "tools/tools.h"
#include "parallel/parallel.h"
#include "biology/DNASequence.h"
#include "scan/SequenceScanner.h"

#include "filesystem"
#include "parallel/parallel.h"
#include <mutex>
// Helper functor to adapt `FILE*` to `std::istream`
class PopenBuffer : public std::streambuf
{
	public:
		PopenBuffer(FILE* file) : file_(file) {}
	protected:
		int underflow() override
		{
			if (gptr() == egptr())
			{
				size_t n = fread(buffer_, 1, sizeof(buffer_), file_);
				if (n == 0) return traits_type::eof();
				setg(buffer_, buffer_, buffer_ + n);
			}
			return traits_type::to_int_type(*gptr());
		}
	private:
		FILE* file_;
		char buffer_[1024];
};


Sequence::DNA ParseRawInput(std::istream& inputStream,SequenceScanner & scanner)
{
	// DataHolder data;
	// StatefulReader Reader(data);

	// int counter = 0;
	std::string PIPE_LINE;
	bool readNextLine = false;
	Sequence::DNA seq("");
	Record rec;
	Sequence::DNA bestseq("");

	double bestScore = -9999999;
	while (std::getline(inputStream,PIPE_LINE))
	{
		if (PIPE_LINE[0]=='@')
		{
			readNextLine = true;
		}
		else 
		{
			if (readNextLine)
			{
				try
				{
					seq.NewSequence(PIPE_LINE);
					scanner.Scan(seq,rec);
					if (bestScore < rec.Score)
					{
						bestScore = rec.Score;
						bestseq = seq;
						// LOG(INFO) << PIPE_LINE << " " << seq.FileString;
					}
				}
				catch (...)
				{
					return bestseq;
				}
			}
			readNextLine = false;
		}
		
		// std::string chromosome;
        // dnaindex idx;
        // unsigned int k;

        // parseLine(PIPE_LINE, Settings.StreamDelimiter, chromosome, idx, k);

		// Reader.AddLine(chromosome,idx,k);
	}
	// if (counter > maxPerFile)
	// {
	// 	maxPerFile = counter;
	// }
	// globalCounter += counter;
	// LOG(INFO) << counter << " for a total " << globalCounter << " max = " << maxPerFile;
	return bestseq;
	
}

Sequence::DNA ShellExecute(std::string cmd,SequenceScanner & scanner)
{
	// std::string cmd = "cat Data/Aaron.dat";
	LOG(DEBUG) << "Calling popen with command '" << cmd << "'";
	FILE* pipe = popen(cmd.c_str(),"r");
	if (!pipe)
	{
		throw std::runtime_error("Failed to open pipe for command: " + cmd);
	}
	PopenBuffer tmp(pipe);
	std::istream stream(&tmp);
	// auto out = ParseRawInput(stream);
	auto out = ParseRawInput(stream,scanner);
	auto exit = pclose(pipe);
	if (WEXITSTATUS(exit) != 0)
	{
		throw std::runtime_error("Command (" + cmd +") returned a non-zero exit code");
	}
	return out;
}

void scanSeq()
{

	
	auto pwm = getRecursiveFileList(Settings.Input.PFMDirectory,Settings.Input.PFMRegex);
	// pwm.resize(3);
	int Nseq = Settings.System.InputSize;
	int SeqLength = 25;
	auto scanner = SequenceScanner(pwm,Nseq,SeqLength);


	auto fastqFiles = getRecursiveFileList(Settings.Input.ReadDirectory,Settings.Input.ReadRegex);
	
	
	fastqFiles.resize(5);
	
	
	ParallelPool Parallel(Settings.System.ParallelThreads);

	std::atomic<int> globalCount;
	std::mutex lock;
	ProgressBar PB(fastqFiles.size(),"Iterating through files\n");
	Parallel.For(fastqFiles.size(),[&fastqFiles,&scanner,&globalCount,&lock,&PB](int i)
	{
		// int count = 0;
		auto s = ShellExecute("gzcat " + fastqFiles[i].path().string(),scanner);
		// LOG(INFO) << fastqFiles[i].path().filename() << " best score = " << s.FileString;
		++globalCount;
		lock.lock();
		PB.Update(globalCount);
		lock.unlock();
		// globalCount += s;
		// LOG(INFO) << i << " " << fastqFiles[i].path();
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
	LOG(INFO) << "here!";
	//Load PWMs from memory

	//Load sequences from memory & pre-analyse are we assuming from file and not piped in? Pipes may be more efficient than saving to disk if preprocessing

	//Then run the analysis

	return 0;
}