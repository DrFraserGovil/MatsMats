#include "recursiveFileSearch.h"
#include <regex>
std::vector<std::filesystem::directory_entry> getRecursiveFileList(const std::string & rootPath,const std::string & regexp)
{
	std::regex pattern(regexp);
	std::vector<std::filesystem::directory_entry> paths;
	for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(rootPath))
	{
		if (dirEntry.is_regular_file())
        {
			if (std::regex_match(dirEntry.path().string(),pattern))
			{
				paths.push_back(dirEntry);
			}
		}
		// std::string p = dirEntry.path();
		// if (p.find("fastq") != std::string::npos)
		// {
		// 	LOG(INFO) << p;
		// }
	}
	return paths;
}