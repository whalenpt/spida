
#include <fstream>
#include <vector>
#include <string>
#include <ios>
#include <filesystem>
#include "pwutils/read/readfile.h"
#include "pwutils/read/dat.hpp"
#include "pwutils/read/json.hpp"
#include "pwutils/pwstrings.h"
#include "pwutils/pwmath.hpp"

namespace pw{

    FileSignature fileSignature(const std::filesystem::path& path){
        std::string file_signature = pw::stringLowerCase(path.extension().string());
        if(file_signature == ".json")
            return FileSignature::JSON;
        else if(file_signature == ".dat")
            return FileSignature::DAT;
        return deduceFileSignature(path);
    }

    FileSignature deduceFileSignature(const std::filesystem::path& path)
    {
        std::ifstream stream{path};
        std::string line;
        while(std::getline(stream,line)){
            line = pw::eatWhiteSpace(line);
            if(line.empty() || line.front() == '#')
                continue;
            else if(line == "{"){
                // JSON? check to be sure
                return checkJSONSignature(stream,line);
            } else if(pw::lineIsIntegers(line) || pw::lineIsDoubles(line)){
                // DAT? check to be sure
                return checkDatSignature(stream,line);
            }
        }
        return FileSignature::UNKNOWN;
    }

    FileSignature checkDatSignature(std::ifstream& stream,std::string& str){
        std::string line = pw::eatWhiteSpace(str);
        std::vector<std::string> line_data = pw::parseString(line,' ');
        if(!pw::rowIsDoubles(line_data))
            return FileSignature::UNKNOWN;

        // Check for next line data
        if(line_data.size() == 2 && pw::rowIsIntegers(line_data)){
            int n1 = std::stoi(line_data[0]);
            int n2 = std::stoi(line_data[1]);
            getDatLineData(stream,line_data);
            if(line_data.size() == 1){ // Check that data is XYZ DAT format
                int count = 0;
                while(count < (n1+n2) && std::getline(stream,line))
                    count++;
                if(stream.eof())
                    return pw::FileSignature::UNKNOWN;
                line_data = pw::parseString(line,' ');
                if(!pw::rowIsDoubles(line_data) || line_data.size() != n1)
                    return pw::FileSignature::UNKNOWN;
                return pw::FileSignature::DAT;
            } else if(line_data.size() == 2) // Assume data is XY format (two-columns)
                return pw::FileSignature::DAT;
        } else if(line_data.size() == 2 && pw::rowIsDoubles(line_data)){
            getDatLineData(stream,line_data);
            if(line_data.size() == 2 && pw::rowIsDoubles(line_data))
                return pw::FileSignature::DAT;
            else 
                return pw::FileSignature::UNKNOWN;
        } else if(line_data.size() == 3 && pw::rowIsDoubles(line_data)){
            getDatLineData(stream,line_data);
            if(line_data.size() == 3 && pw::rowIsDoubles(line_data))
                return pw::FileSignature::DAT;
            else 
                return pw::FileSignature::UNKNOWN;
        } 
        return pw::FileSignature::UNKNOWN;
    }

    void getDatLineData(std::ifstream& stream,std::vector<std::string>& line_data)
    {
        std::string line;
        while(std::getline(stream,line)){
           line = pw::eatWhiteSpace(line);
           // Ignore empty lines and comment lines
           if(line.empty() || line.front() == '#')
               continue;
           else{
               // BINGO - split by empty spaces
               line_data = pw::parseString(line,' ');
               return;
           }
        }
        throw std::domain_error("No DAT data found in stream");
    }

    FileSignature checkJSONSignature(std::ifstream& stream,std::string& line){
        if(std::getline(stream,line)){
            line  = pw::eatWhiteSpace(line);
            // Check the first line to see if its in JSON form
            if(!checkJSONline(line))
                return FileSignature::UNKNOWN;
            if(&line.back() == std::string(",")){
                // Check second line
                if(std::getline(stream,line))
                    if(checkJSONline(line))
                        return FileSignature::JSON;
            } else
                return FileSignature::JSON;
        }
        return FileSignature::UNKNOWN;
    }

    bool checkJSONline(std::string& line){
        line  = pw::eatWhiteSpace(line);
        std::vector<std::string> colon_data = pw::parseString(line,':');
        if(colon_data.size() != 2)
            // Expect one colon per JSON item
            return false;
        if(colon_data[0].at(0) != '"')
            // Expect Parentheses to start JSON label
            return false;
        return true;
    }

    DataSignature dataSignature(const std::filesystem::path& path,FileSignature file_signature)
    {
        if(file_signature == FileSignature::DAT)
            return dat::dataSignature(path);
        else if(file_signature == FileSignature::JSON)
            return json::dataSignature(path);
        else
            return DataSignature::UNKNOWN;
    }


    OperatorSignature operatorSignature(const std::filesystem::path& path,FileSignature file_signature)
    {
        if(file_signature == FileSignature::DAT)
            return dat::operatorSignature(path);
        else if(file_signature == FileSignature::JSON)
            return json::operatorSignature(path);
        else
            return OperatorSignature::NONE;
    }



}






