
#include "pwutils/report/basereport.h"
#include <string>
#include <map>

namespace pw{

void ReportBase::setItem(const std::string& key,double val) {
    std::string strVal = std::to_string(val);
    m_metadata_map[key] = strVal;
}

void ReportBase::setItem(const std::string& key,const std::string& nm) {
  m_metadata_map[key] = nm;
}

void ReportBase::removeItem(const std::string& key) {
  m_metadata_map.erase(key);
}


}

