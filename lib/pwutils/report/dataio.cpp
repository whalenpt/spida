
#include "pwutils/report/dataio.hpp"
#include "pwutils/report/reporthelper.h"
#include "pwutils/pwexcept.h"
#include <filesystem>
#include <stdexcept>
#include <fstream>

namespace pw{

DataIO::DataIO(const fs::path& dir_path) : m_dirpath(dir_path)
{
    createDirectory(m_dirpath,false);
}

void DataIO::setDirectoryPath(const fs::path& dirpath)
{
    if(fs::is_directory(dirpath)){
       m_dirpath = dirpath;
    } else{
        std::string str("DataIO::DataIO(std::filesystem::path& path) \
             specified directory path IS NOT a directory path ");
        throw pw::Exception(str);
    }
}

void DataIO::checkExistence(const fs::path& filepath) const{
    if(!fs::exists(filepath))
        throw std::runtime_error("Filepath of " + filepath.string() + " does not exist");
}

void DataIO::checkOpens(const fs::path& filepath) const{
    std::ifstream fin;
    fin.open(filepath);
    if(!fin.is_open())
        throw std::runtime_error("Could not open file: " + filepath.string());
    fin.close();
}





}







