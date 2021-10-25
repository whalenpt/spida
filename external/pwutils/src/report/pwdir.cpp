

#include <cstring>
#include <fstream>
#include <string>
#include <iostream>
#include <filesystem>
#include <sstream>
#include "pwutils/pwstrings.h"
#include "pwutils/pwexcept.h"
#include "pwutils/report/pwdir.h"

namespace pw{

DirAux::DirAux()
{
	current_dir = std::filesystem::current_path();
}

std::filesystem::path DirAux::stringToPath(std::string dir_name)
{
	namespace fs = std::filesystem;
	std::string local_path = pw::eatWhiteSpace(dir_name);
	local_path = pw::eatWhiteSpace(local_path," \t","_");
	if(local_path.empty())
			throw pw::Exception("Invalid string provided to createDirectory function.");
	fs::path local_dir = fs::path(local_path);
	fs::path full_path = current_dir / local_path;
	return full_path;
}


bool DirAux::checkDirectoryExists(std::filesystem::path dir_path)
{
	if(std::filesystem::exists(dir_path))
			return true;
	return false;
}

bool DirAux::checkDirectoryExists(std::string dir_name)
{
	std::filesystem::path dir_path = stringToPath(dir_name);
	return checkDirectoryExists(dir_path);
}

void DirAux::createDirectory(std::filesystem::path dir_path,bool overwrite)
{
	namespace fs = std::filesystem;
	if(fs::exists(dir_path) && overwrite){
        fs::remove_all(dir_path);
	}
	if(!fs::exists(dir_path)){
        fs::create_directory(dir_path);
	}
}

void DirAux::createDirectory(std::string dir_name,bool overwrite)
{
	std::filesystem::path dir_path = stringToPath(dir_name);
	createDirectory(dir_path,overwrite);
}

std::filesystem::path DirAux::addTargetToDir(std::string target_name,std::filesystem::path dir_path)
{
	namespace fs = std::filesystem;
	if(!checkDirectoryExists(dir_path))
			throw pw::Exception("Invalid directory path provided in DirAux::addTargetToDir");
	fs::path target_path = dir_path / fs::path(target_name);
	return target_path;
}

std::filesystem::path DirAux::addTargetToDir(std::string target_name,std::string dir_name)
{
	std::filesystem::path dir_path = stringToPath(dir_name);
	return addTargetToDir(target_name,dir_path);
}

}




