//readfile.h
#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include "pwutils/pwdefs.h"

namespace pw{
    FileSignature fileSignature(const std::filesystem::path& path);
    FileSignature deduceFileSignature(const std::filesystem::path& path);
    FileSignature checkJSONSignature(std::ifstream& fin,std::string& line);
    FileSignature checkDatSignature(std::ifstream& fin,std::string& line);
    bool checkJSONline(std::string& line);
    void getDatLineData(std::ifstream& infile,std::vector<std::string>& line_data);

    DataSignature dataSignature(const std::filesystem::path& path,FileSignature file_signature);
    OperatorSignature operatorSignature(const std::filesystem::path& path,FileSignature file_signature);
}




