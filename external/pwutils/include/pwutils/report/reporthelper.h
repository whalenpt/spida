
#ifndef REPORT_HELPER_H_
#define REPORT_HELPER_H_

#include <vector>
#include <map>
#include <string>
#include <filesystem>
#include "pwutils/pwconstants.h"

namespace pw{


const int REPORT_PRECISION = 16;
const int REPORT_PADING = 8;
const std::string DEFAULT_REPORTOUT_DIR("pwout");
using metadataMap = std::map<std::string,std::string>;
using metadataPair = std::pair<std::string,std::string>;

namespace fs = std::filesystem;
fs::path createDirectory(const std::string& dir_name, bool overwrite = true);
void createDirectory(const std::filesystem::path& dir_path,bool overwrite = true);
void clearDirectory(const std::filesystem::path& dir_path);

fs::path filePath(const fs::path& dir_path,\
        const std::string& nm,int repNum,const std::string& extension); 
fs::path filePath(const fs::path& dir_path,\
        const std::string& nm,const std::string& extension);

void getPhaseLimits(std::vector<double>& in,int& stindx,int& endindx,int sz);
bool CheckSignChange(std::vector<double>& in,int indx1,int indx2);
void AdjustPhase(std::vector<double>& in,int sz);
int ComputeOutN_2D(int nD,int stride);

}

#endif

