#pragma once


#include "Log.h"


template<char Symbol = '#',unsigned int MaxMarks =16>
class ProgressBar
{
	private:
int Count;
		std::string Prefix;
		std::string Suffix;
		int prevMarks;
		bool IsActive;
		double marksPerProgress;

		
	public:
		ProgressBar(int maxCount,std::string prefix="", std::string suffix="") : Count(maxCount), Prefix(prefix), Suffix(suffix)
		{
			IsActive = GlobalLog::Config.TerminalOutput;
			prevMarks = -1;
			marksPerProgress = MaxMarks * 1.0/(maxCount-1); 
			
			if (!IsActive)
			{
				LOG(INFO) << Prefix << "[~~Progress Bar disabled~~]" << Suffix;
			}
		}

		void Update(int progress)
		{
			if (IsActive)
			{
				int newMarks = std::min((uint)round(marksPerProgress * progress),MaxMarks);
				
				if (newMarks > prevMarks || prevMarks < 0)
				{
					if (prevMarks >= 0)
					{
						LOG(INFO).ErasePrevious();
					}
					LOG(INFO) << Prefix << "[" << std::string(newMarks,Symbol) << std::string(MaxMarks - newMarks,' ') << "]" << Suffix;
					prevMarks = newMarks;
				}

			}
		}
};


namespace Progress
{
	template<class T>
	void inline For(int start, int end,T callable,std::string pbPrefix= "", std::string pbSuffix="")
	{
		ProgressBar pb(end-start,pbPrefix,pbSuffix);
		for (int i =start; i < end; ++i)
		{
			pb.Update(i);
			callable(i);
		}
	}
}