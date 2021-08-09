
#ifndef DATAIO_HPP_
#define DATAIO_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include "pwutils/report/reporthelper.h"

namespace pw{

namespace fs = std::filesystem;

class DataIO{
    public:
        DataIO(const fs::path& dirpath = fs::current_path());
        void setDirectoryPath(const fs::path& dirpath);
        template<typename T>
        void writeFile(const fs::path& fname,const std::vector<T>& x) const;
        template<typename T1,typename T2>
        void writeFile(const fs::path& fname,const std::vector<T1>& x,const std::vector<T2>& y) const;

        // read one col data
        template<typename T>
        void readFile(const fs::path& fname,std::vector<T>&x) const;
        // read two col data
        template<typename T1,typename T2>
        void readFile(const fs::path& fname,std::vector<T1>&x,std::vector<T2>&y) const;
        // read three col data
        template<typename T1,typename T2,typename T3>
        void readFile(const fs::path& fname,std::vector<T1>&x,std::vector<T2>&y,std::vector<T3>&z) const;

        template<typename T>
        void appendFile(const fs::path& fname,const T& val) const;
        template<typename T1,typename T2>
        void appendFile(const fs::path& fname,const T1& val1,const T2& val2) const;

        void clearDirectory() { pw::clearDirectory(m_dirpath); }
    private:
        fs::path m_dirpath;
};

template<typename T>
void DataIO::writeFile(const fs::path& fname,const std::vector<T>& x) const
{
    fs::path file_path = m_dirpath / fname;
    std::ofstream fout(file_path);
    for(auto val : x)
        fout << std::scientific << std::setprecision(pw::REPORT_PRECISION) 
             << std::setw(pw::REPORT_PRECISION + pw::REPORT_PADING) << val << std::endl;
    fout.close();
}

template<typename T1,typename T2>
void DataIO::writeFile(const fs::path& fname,const std::vector<T1>& x,const std::vector<T2>& y) const
{
    assert(x.size() == y.size());
    fs::path file_path = m_dirpath / fname;
    std::ofstream fout(file_path);
    for(auto j = 0; j < x.size(); j++){
        fout << std::scientific << std::setprecision(pw::REPORT_PRECISION) 
            << std::setw(pw::REPORT_PRECISION + pw::REPORT_PADING) << x[j] 
            << std::scientific << std::setprecision(pw::REPORT_PRECISION) 
            << std::setw(pw::REPORT_PRECISION + pw::REPORT_PADING) << y[j] << std::endl; 
    }
    fout.close();
}

template<typename T>
void DataIO::appendFile(const fs::path& fname,const T& val) const
{
    fs::path file_path = m_dirpath / fname;
    std::ofstream fout(file_path,std::ofstream::out | std::ofstream::app);
    fout << std::scientific << std::setprecision(pw::REPORT_PRECISION) 
         << std::setw(pw::REPORT_PRECISION + pw::REPORT_PADING) << val << std::endl;
    fout.close();
}

template<typename T1,typename T2>
void DataIO::appendFile(const fs::path& fname,const T1& val1,const T2& val2) const
{
    fs::path file_path = m_dirpath / fname;
    std::ofstream fout(file_path,std::ofstream::out | std::ofstream::app);
    fout << std::scientific << std::setprecision(pw::REPORT_PRECISION) 
         << std::setw(pw::REPORT_PRECISION + pw::REPORT_PADING) << val1 
         << std::scientific << std::setprecision(pw::REPORT_PRECISION) 
         << std::setw(pw::REPORT_PRECISION + pw::REPORT_PADING) << val2 << std::endl; 
    fout.close();
}

template<typename T>
void DataIO::readFile(const fs::path& fname,std::vector<T>&x) const
{
    fs::path file_path = m_dirpath / fname;
    std::ifstream fin(file_path);
    T val;
    while(fin >> val)
        x.push_back(val);
    fin.close();
}

template<typename T1,typename T2>
void DataIO::readFile(const fs::path& fname,std::vector<T1>&x,std::vector<T2>&y) const
{
    fs::path file_path = m_dirpath / fname;
    std::ifstream fin(file_path);
    T1 val1; T2 val2;
    while(fin >> val1 >> val2){
        x.push_back(val1);
        y.push_back(val2);
    }
    fin.close();
}

template<typename T1,typename T2,typename T3>
void DataIO::readFile(const fs::path& fname,std::vector<T1>&x,std::vector<T2>&y,std::vector<T3>&z) const
{
    fs::path file_path = m_dirpath / fname;
    std::ifstream fin(file_path);
    T1 val1; T2 val2; T3 val3;
    while(fin >> val1 >> val2 >> val3){
        x.push_back(val1);
        y.push_back(val2);
        z.push_back(val3);
    }
    fin.close();
}








}

#endif

