
#include "pwutils/basereport.h"
#include "pwutils/reporthelper.h"
#include <string>
#include <map>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace pw{

void ReportBase::report(std::ofstream& os) const
{
    pw::createDirectory(m_dirpath, false);
    os.open(ReportBase::path());
    if(!os.is_open())
        throw std::runtime_error("Failed to open data file for output stream");
    reportImplement(os);
    os.close();
}

void ReportBase::report(std::ofstream& os, unsigned rep_num) const
{
    pw::createDirectory(m_dirpath, false);
    os.open(ReportBase::path(rep_num));
    if(!os.is_open())
        throw std::runtime_error("Failed to open data file for output stream");
    reportImplement(os);
    os.close();
}

std::ofstream& operator<<(std::ofstream& os, const ReportBase& def) {
    def.report(os);
    return os;
}

std::ofstream& operator<<(std::ofstream& os, const ReportBase* def) {
    def->report(os);
    return os;
}

std::ofstream& operator<<(std::ofstream& os, const std::unique_ptr<ReportBase>& def) {
    def->report(os);
    return os;
}

void ReportBase::setItem(const std::string& key, double val) {
    std::ostringstream oss;
    oss << std::setprecision(12) << val;
    m_metadata_map[key] = oss.str();
}

void ReportBase::setItem(const std::string& key, float val) {
    std::ostringstream oss;
    oss << std::setprecision(6) << val;
    m_metadata_map[key] = oss.str();
}

void ReportBase::setItem(const std::string& key, int val) {
    std::ostringstream oss;
    oss << val;
    m_metadata_map[key] = oss.str();
}

void ReportBase::setItem(const std::string& key, const std::string& nm) { m_metadata_map[key] = nm; }
void ReportBase::removeItem(const std::string& key) { m_metadata_map.erase(key); }

}
