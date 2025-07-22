#pragma once

#include <vector>
#include <string>
#include <filesystem>
std::vector<std::filesystem::directory_entry> getRecursiveFileList(const std::string &  rootPath,const std::string & regexp = ".");