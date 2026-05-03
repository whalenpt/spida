
#include "pwutils/report/basereport.h"
#include "pwutils/report/reporthelper.h"
#include <string>
#include <map>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace pw{

void ReportBase::report(std::ofstream& os) const
{
    pw::createDirectory(m_dirpath,false);
    os.open(ReportBase::path());
    if(!os.is_open())
        throw std::runtime_error("Failed to open data file for output stream");
    reportImplement(os);
    os.close();
}

void ReportBase::report(std::ofstream& os,unsigned rep_num) const
{
    pw::createDirectory(m_dirpath,false);
    os.open(ReportBase::path(rep_num));
    if(!os.is_open())
        throw std::runtime_error("Failed to open data file for output stream");
    reportImplement(os);
    os.close();
}


std::ofstream& operator<<(std::ofstream& os,const ReportBase& def){
    def.report(os);
    return os;
}

std::ofstream& operator<<(std::ofstream& os,const ReportBase* def){ 
    def->report(os);
    return os;
}

std::ofstream& operator<<(std::ofstream& os,const std::unique_ptr<ReportBase> def)
{ 
    def->report(os);
    return os;
}

void ReportBase::setItem(const std::string& key,double val) {
    std::ostringstream oss;
    oss << std::setprecision(12) << val;
    m_metadata_map[key] = oss.str();
}

void ReportBase::setItem(const std::string& key,float val) {
    std::ostringstream oss;
    oss << std::setprecision(6) << val;
    m_metadata_map[key] = oss.str();
}

void ReportBase::setItem(const std::string& key,int val) {
    std::ostringstream oss;
    oss << std::setprecision(12) << val;
    m_metadata_map[key] = oss.str();
}


void ReportBase::setItem(const std::string& key,const std::string& nm) { m_metadata_map[key] = nm; }
void ReportBase::removeItem(const std::string& key) { m_metadata_map.erase(key); }

void ReportBase::setFileSignature(FileSignature file_sig) {
    m_metadata_map["FileSignature"] = std::to_string(static_cast<int>(file_sig));
}

void ReportBase::setDataSignature(DataSignature data_sig){
    m_metadata_map["DataSignature"] = std::to_string(static_cast<int>(data_sig));
}
void ReportBase::setOperatorSignature(OperatorSignature op_sig) {
    m_metadata_map["OperatorSignature"] = std::to_string(static_cast<int>(op_sig));
}

FileSignature ReportBase::fileSignature() const
{
    const auto it = m_metadata_map.find("FileSignature");
    if(it != m_metadata_map.end())
        return static_cast<FileSignature>(std::stoi(it->second));
    return FileSignature::UNKNOWN;
}

DataSignature ReportBase::dataSignature() const
{
    const auto it = m_metadata_map.find("DataSignature");
    if(it != m_metadata_map.end())
        return static_cast<DataSignature>(std::stoi(it->second));
    return DataSignature::UNKNOWN;
}

OperatorSignature ReportBase::operatorSignature() const
{
    const auto it = m_metadata_map.find("OperatorSignature");
    if(it != m_metadata_map.end())
        return static_cast<OperatorSignature>(std::stoi(it->second));
    return OperatorSignature::NONE;
}


}

