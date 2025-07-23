#pragma once
#include "filesystem"
#include <streambuf>
#include "SequenceScanner.h"
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


void inline parseLine(std::string_view fileLine,SequenceScanner & scanner, Sequence::DNA & dna, Record & record, bool & nextLineFlag,std::string & gatheredID,std::ofstream & file)
{
	if (fileLine[0]=='@')
	{
		nextLineFlag = true; 
		auto firstSpace = std::find(fileLine.begin(),fileLine.end(),' ');
		gatheredID = std::string(fileLine.begin()+1,firstSpace);
	}
	else 
	{
		if (nextLineFlag)
		{	
			dna.NewSequence(fileLine);
			if (dna.AlphabetContained)
			{
				scanner.Scan(dna,record);
				file << gatheredID << " " << dna.FileString << "\n";
			}
		}
		nextLineFlag = false;
	}
}



//this calls to an external tool (gzcat), which unzips the file and then spits it out for us to catch
void gzfastQScan(const std::string & filename,SequenceScanner & scanner,std::ofstream & file)
{
	std::string cmd = "gzcat " + filename;
	LOG(DEBUG) << "Calling popen with command '" << cmd << "'";
	
	FILE* pipe = popen(cmd.c_str(),"r");
	if (!pipe)
	{
		throw std::runtime_error("Failed to open pipe for command: " + cmd);
	}
	
	PopenBuffer tmp(pipe);
	std::istream stream(&tmp);
	
	std::string PIPE_LINE;
	
	Sequence::DNA seq("");
	Record rec;
	std::string previousLine;
	bool readNextLine = false;
	while (std::getline(stream,PIPE_LINE))
	{
		parseLine(PIPE_LINE,scanner,seq,rec,readNextLine,previousLine,file);
	}

	auto exit = pclose(pipe);
	if (WEXITSTATUS(exit) != 0)
	{
		throw std::runtime_error("Command (" + cmd +") returned a non-zero exit code");
	}
}

void fastqScan(const std::string & filename,SequenceScanner & scanner,std::ofstream & file)
{
	Sequence::DNA seq("");
	Record rec;
	bool readNextLine = false;
	std::string previousLine;
	forLineIn(filename,[&](auto line){
		parseLine(line,scanner,seq,rec,readNextLine,previousLine,file);
	});
}