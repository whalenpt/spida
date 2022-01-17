// dat.hpp
#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <complex>
#include "pwutils/pwdefs.h"

namespace dat{
    using pw::dcmplx;
    pw::DataSignature dataSignature(const std::filesystem::path& path);
    pw::DataSignature deduceDataSignature(std::ifstream& fin);

    pw::OperatorSignature operatorSignature(const std::filesystem::path& path);
    pw::metadataMap getHeaderContent(std::ifstream& iss);
    void getLineOfData(std::ifstream& fin,std::vector<std::string>& line_data);

    template<typename T1,typename T2>
    pw::metadataMap readXY(const std::filesystem::path& path,std::vector<T1>& x,
            std::vector<T2>& y)
    {
        std::ifstream infile{path};
        pw::metadataMap metadata = getHeaderContent(infile);
        while(!infile.eof()){
            T1 a; 
            T2 b;
            infile >> a >> b;
            x.push_back(a);
            y.push_back(b);
        }
        return metadata;
    } 
 

    template<typename T1,typename T2>
    pw::metadataMap readXCVY(const std::filesystem::path& path,std::vector<T1>& x,
            std::vector<std::complex<T2>>& y)
    {
        std::ifstream infile{path};
        pw::metadataMap metadata = getHeaderContent(infile);
        while(!infile.eof()){
            T1 a;
            T2 b,c;
            infile >> a >> b >> c;
            x.push_back(a);
            y.push_back(std::complex<T2>(b,c));
        }
        return metadata;
    } 

    template<typename T1,typename T2>
    void readXY3D(std::ifstream& fin,std::vector<T1>& x,std::vector<T2>& y)
    {
        unsigned int nx,ny;
        fin >> nx >> ny;
        x.resize(nx);
        y.resize(ny);
        for(auto i = 0; i < nx; i++)
            fin >> x[i];
        for(auto i = 0; i < ny; i++)
            fin >> y[i];
    }
     
    template<typename T1,typename T2,typename T3>
    pw::metadataMap readXYZ(const std::filesystem::path& path,std::vector<T1>& x,\
            std::vector<T2>& y,std::vector<T3>& z)
    {
        std::ifstream infile{path};
        pw::metadataMap metadata = getHeaderContent(infile);
        readXY3D<T1,T2>(infile,x,y);
        z.reserve(x.size()*y.size());
        std::string line;
        while(std::getline(infile,line)){
            double val;
            std::stringstream ss(line);
            while(ss >> val)
                z.push_back(val);
        }
        return metadata;
    } 
     



    template<typename T1,typename T2,typename T3>
    pw::metadataMap readXYCVZ(const std::filesystem::path& path,std::vector<T1>& x,\
            std::vector<T2>& y,std::vector<std::complex<T3>>& z)
    {
        std::ifstream infile{path};
        pw::metadataMap metadata = getHeaderContent(infile);
        readXY3D<T1,T2>(infile,x,y);
        z.resize(x.size()*y.size());
        std::string line;
        while(std::getline(infile,line)){
            T3 real_val, imag_val;
            std::stringstream ss(line);
            while(ss >> real_val >> imag_val)
                z.push_back(dcmplx(real_val,imag_val));
        }
        return metadata;
    }


}



