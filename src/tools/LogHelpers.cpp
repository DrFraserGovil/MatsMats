#include "LogHelpers.h"
#include "Log.h"

bool isTerminal() {
    return isatty(fileno(stdout));
}

namespace GlobalLog
{
	
	ConfigObject::ConfigObject()
	{
		SetLevel(INFO);
		ShowHeaders = true;
		AppendNewline = true;
		TerminalOutput = isTerminal();
	}

	void ConfigObject::SetLevel(int level)
	{
		Level = LogLevelConvert(level);
	}
	

	void ConfigObject::Initialise(int level, bool header)
	{
		AppendNewline =true;
		SetLevel(level);
		
		ShowHeaders = header;
		LOG(DEBUG) << "Logging system initialised";
		
	}


	//See the LogHelpers.h file for the definitions of these global variables
	ConfigObject Config; 
	std::vector<size_t> PreviousLines = std::vector<size_t>(LogLevel::MAXLEVEL,0);
	std::mutex StreamMutex;
}