
#include "pwutils/report/dataio.hpp"
#include "pwutils/report/reporthelper.h"
#include "pwutils/pwexcept.h"
#include <filesystem>
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




}







