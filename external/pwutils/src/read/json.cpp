
#include "pwutils/read/json.hpp"
#include "pwutils/pwstrings.h"
#include "pwutils/pwmath.hpp"
#include <fstream>
#include <string>
#include <map>
#include <stdexcept>
#include <json11/json11.hpp>

namespace json{


pw::metadataMap getMetaData(const std::filesystem::path& path)
{
    std::ifstream stream{path};
    return getMetaData(stream);
}

pw::metadataMap getMetaData(std::ifstream& iss){

    size_t length = 16;
    char* buffer = new char[length+1];
    std::string metastring;
    bool search = true;
    // read until first array is found
    while(search && iss.read(buffer,length)){
        buffer[length] = '\0';
        std::string temp(buffer);
        auto found = temp.find_first_of('[');
        if(found != std::string::npos){
            metastring += temp.substr(0,found);
            search = false;
        }
        else
            metastring += temp;
    }
    // iss reached end
    if(!iss){
        std::string temp(buffer);
        auto found = temp.find_first_of('[');
        if(found != std::string::npos)
            metastring += temp.substr(0,found);
        else
            throw std::domain_error("Failed to find any vector data in the JSON file");
    }
    auto index  = metastring.find_last_of(',');
    metastring = metastring.substr(0,index);
    metastring += '\n';
    metastring += '}';

    delete [] buffer;
    
    std::string err_str;
    const auto meta_json = json11::Json::parse(metastring,err_str);
    const auto& json_map = meta_json.object_items();
    pw::metadataMap metamap;
    for(const auto& item : json_map)
        metamap.insert(std::pair(item.first,item.second.string_value()));
    return metamap;
}

pw::DataSignature dataSignature(const std::filesystem::path& path) {
    std::ifstream stream{path};
    pw::metadataMap meta_map = getMetaData(stream);
    if(meta_map.find("DataSignature") != meta_map.end())
        return static_cast<pw::DataSignature>(std::stoi(meta_map["DataSignature"]));
    else
        return deduceDataSignature(stream);
}

void readJSONstring(const std::filesystem::path& path,std::string& str)
{
    std::ifstream stream{path};
    stream.clear();
    stream.seekg(0,std::ios::end);
    str.resize(stream.tellg());
    stream.seekg(0,std::ios::beg);
    stream.read(&str[0],str.size());
    stream.close();
}

void readJsonObject(const std::filesystem::path& path,json11::Json& json_obj)
{
    std::string buffer;
    readJSONstring(path,buffer);
    std::string err_str;
    json_obj = json11::Json::parse(buffer,err_str);
}

pw::DataSignature deduceDataSignature(std::ifstream& fin)
{
    fin.clear();
    fin.seekg(0,std::ios::end);
    std::string buffer;
    buffer.resize(fin.tellg());
    fin.seekg(0,std::ios::beg);
    fin.read(&buffer[0],buffer.size());
    fin.close();

    std::string err_str;
    const auto json_obj = json11::Json::parse(buffer,err_str);
    if(!err_str.empty())
        throw std::runtime_error(err_str);
    const auto& json_map = json_obj.object_items();

    
    bool xfound = false;
    bool yfound = false;
    bool zfound = false;
    bool data_complex = false;

    for(const auto& item : json_map){
        std::string key = item.first;
//        std::cout << key << std::endl;
        if(key == pw::XLABEL && item.second.is_array())
            xfound = true;
        else if(key == pw::YLABEL && item.second.is_array())
            yfound = true;
        else if(key == pw::ZLABEL && item.second.is_array())
            zfound = true;
        else if(key == "dtype" && item.second.is_string()){
            if(item.second.string_value() == "complex")
                data_complex = true;
        }
    }

    if(xfound && yfound && !zfound){
        if(!data_complex)
            return pw::DataSignature::XY;
        return pw::DataSignature::XCVY;
    } else if(xfound && yfound && zfound){
        if(!data_complex)
            return pw::DataSignature::XYZ;
        return pw::DataSignature::XYCVZ;
    }
    return pw::DataSignature::UNKNOWN;
}

pw::OperatorSignature operatorSignature(const std::filesystem::path& path)
{
    std::ifstream stream{path};
    pw::metadataMap meta_map = getMetaData(stream);
    if(meta_map.find("OperatorSignature") != meta_map.end())
        return static_cast<pw::OperatorSignature>(std::stoi(meta_map["OperatorSignature"]));
    else
        return pw::OperatorSignature::NONE;
}

void dataNotFound(const std::string& id)
{
    const std::string str("Error in readVecData: id "+ id + " was not found in the json object"); 
            throw std::domain_error(str);
}

void readVecData(const json11::Json& json_obj,std::vector<double>& vec,const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        vec.resize(json_arr.size());
        for(auto i = 0; i < json_arr.size(); i++)
            vec[i] = json_arr[i].number_value();
    }else
        dataNotFound(id);
}

void readVecData(const json11::Json& json_obj,std::vector<float>& vec,const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        vec.resize(json_arr.size());
        for(auto i = 0; i < json_arr.size(); i++)
            vec[i] = json_arr[i].number_value();
    } else
        dataNotFound(id);
}

void readVecData(const json11::Json& json_obj,std::vector<int>& vec,const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        vec.resize(json_arr.size());
        for(auto i = 0; i < json_arr.size(); i++)
            vec[i] = json_arr[i].int_value();
    } else
        dataNotFound(id);
}

void readVecData(const json11::Json& json_obj,std::vector<std::string>& vec,const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        vec.resize(json_arr.size());
        for(auto i = 0; i < json_arr.size(); i++)
            vec[i] = json_arr[i].string_value();
    } else
        dataNotFound(id);
}

void readVecData(const json11::Json& json_obj,std::vector<std::complex<double>>& vec,\
        const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        if(json_arr.size() % 2 != 0)
            throw std::domain_error("readComplexVecData requires an even number of components in the json array");
        vec.resize(json_arr.size()/2);
        for(auto i = 0; i < json_arr.size()/2;i++)
            vec[i] = std::complex<double>(json_arr[2*i].number_value(),json_arr[2*i+1].number_value());
    } else
        dataNotFound(id);
}

void readVecData(const json11::Json& json_obj,std::vector<std::complex<float>>& vec,\
        const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        if(json_arr.size() % 2 != 0)
            throw std::domain_error("readComplexVecData requires an even number of components in the json array");
        vec.resize(json_arr.size()/2);
        for(auto i = 0; i < json_arr.size()/2;i++)
            vec[i] = std::complex<float>(static_cast<float>(json_arr[2*i].number_value()),\
                    static_cast<float>(json_arr[2*i+1].number_value()));
    } else
        dataNotFound(id);
}

void readVecData(const json11::Json& json_obj,std::vector<std::complex<int>>& vec,\
        const std::string& id)
{
    if(json_obj[id].is_array()){
        const json11::Json::array& json_arr = json_obj[id].array_items();
        if(json_arr.size() % 2 != 0)
            throw std::domain_error("readComplexVecData requires an even number of components in the json array");
        vec.resize(json_arr.size()/2);
        for(auto i = 0; i < json_arr.size()/2;i++)
            vec[i] = std::complex<int>(json_arr[2*i].int_value(),json_arr[2*i+1].int_value());
    
    } else
        dataNotFound(id);
}



}











