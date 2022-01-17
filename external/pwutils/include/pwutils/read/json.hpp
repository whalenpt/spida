// json.hpp
#pragma once 

#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <json11/json11.hpp>
#include "pwutils/pwdefs.h"

namespace json11{
    class Json;
}

namespace json{
    using pw::dcmplx;
    pw::metadataMap getMetaData(std::ifstream& iss);
    pw::metadataMap getMetaData(const std::filesystem::path& path);
    pw::DataSignature dataSignature(const std::filesystem::path& path);
    pw::DataSignature deduceDataSignature(std::ifstream& fin);
    void readJSONstring(const std::filesystem::path& path,std::string& str);
    void readJsonObject(const std::filesystem::path& path,json11::Json& json_obj);

    void readVecData(const json11::Json& json_obj,std::vector<double>& vec,const std::string& id);
    void readVecData(const json11::Json& json_obj,std::vector<float>& vec,const std::string& id);
    void readVecData(const json11::Json& json_obj,std::vector<int>& vec,const std::string& id);
    void readVecData(const json11::Json& json_obj,std::vector<std::string>& vec,const std::string& id);
    void readVecData(const json11::Json& json_obj,std::vector<std::complex<double>>& vec,\
            const std::string& id);
    void readVecData(const json11::Json& json_obj,std::vector<std::complex<float>>& vec,\
            const std::string& id);
    void readVecData(const json11::Json& json_obj,std::vector<std::complex<int>>& vec,\
            const std::string& id);
    void dataNotFound(const std::string& id);

    pw::OperatorSignature operatorSignature(const std::filesystem::path& path);

    template<typename T1,typename T2>
    pw::metadataMap readXY(const std::filesystem::path& path,std::vector<T1>& x,\
            std::vector<T2>& y)
    {
        pw::metadataMap metadata = getMetaData(path);
        json11::Json json_obj;
        readJsonObject(path,json_obj);
        readVecData(json_obj,x,pw::XLABEL);
        readVecData(json_obj,y,pw::YLABEL);
        return metadata;
    } 
     
    template<typename T1,typename T2>
    pw::metadataMap readXCVY(const std::filesystem::path& path,std::vector<T1>& x,\
             std::vector<std::complex<T2>>& y)
    {
        pw::metadataMap metadata = getMetaData(path);
        json11::Json json_obj;
        readJsonObject(path,json_obj);
        readVecData(json_obj,x,pw::XLABEL);
        readVecData(json_obj,y,pw::YLABEL);
        return metadata;
    } 
     

    template<typename T1,typename T2,typename T3>
    pw::metadataMap readXYZ(const std::filesystem::path& path,std::vector<T1>& x,\
            std::vector<T2>& y,std::vector<T3>& z)
    {
        pw::metadataMap metadata = getMetaData(path);
        json11::Json json_obj;
        readJsonObject(path,json_obj);
        readVecData(json_obj,x,pw::XLABEL);
        readVecData(json_obj,y,pw::YLABEL);
        readVecData(json_obj,z,pw::ZLABEL);
        return metadata;
    } 
 
    template<typename T1,typename T2,typename T3>
    pw::metadataMap readXYCVZ(const std::filesystem::path& path,std::vector<double>& x,\
             std::vector<double>& y,std::vector<std::complex<T3>>& z)
    {
        pw::metadataMap metadata = getMetaData(path);
        json11::Json json_obj;
        readJsonObject(path,json_obj);
        readVecData(json_obj,x,pw::XLABEL);
        readVecData(json_obj,y,pw::YLABEL);
        readVecData(json_obj,z,pw::ZLABEL);
        return metadata;
    } 
 
}


