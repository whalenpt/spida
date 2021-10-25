
#include "pwutils/read/dat.hpp"
#include "pwutils/pwdefs.h"
#include "pwutils/pwstrings.h"
#include "pwutils/pwmath.hpp"
#include <string>
#include <map>

namespace dat{

pw::DataSignature dataSignature(const std::filesystem::path& path) {
    std::ifstream stream{path};
    pw::metadataMap meta_map = getHeaderContent(stream);
    if(meta_map.find("DataSignature") != meta_map.end())
        return static_cast<pw::DataSignature>(std::stoi(meta_map["DataSignature"]));
    else
        return dat::deduceDataSignature(stream);
}

pw::DataSignature deduceDataSignature(std::ifstream& fin)
{
    std::vector<std::string> first_line,second_line,third_line;
    getLineOfData(fin,first_line);
    getLineOfData(fin,second_line);
    getLineOfData(fin,third_line);
    if(first_line.empty() || second_line.empty() || third_line.empty())
        return pw::DataSignature::UNKNOWN;
    auto sz1 = first_line.size();
    auto sz2 = second_line.size();
    auto sz3 = third_line.size();
    if((sz1 == sz2) && (sz2 == sz3)){
        // Check that sizes are 2 or 3
        if((sz1 != 2) && (sz1 != 3))
            return pw::DataSignature::UNKNOWN;
        // Check that all the data is doubles (or perhaps floats, but not ints)
        if(!(pw::rowIsDoubles(first_line) && pw::rowIsDoubles(second_line) && pw::rowIsDoubles(third_line)))
            return pw::DataSignature::UNKNOWN;
        // Assume this is two column data throughout
        if(sz1 == 2)
            return pw::DataSignature::XY;
        // Assume this is three column data throughout
        else if(sz1 == 3)
            return pw::DataSignature::XCVY;
    } else{
        std::cout << "Reading dat data (2D?)" << std::endl;
        // Check for XYZ data (first_line has format [nx ny], next nx+ny lines are size 1
        if(sz1 != 2 || sz2 != 1 || sz3 != 1)
            return pw::DataSignature::UNKNOWN;
        int nx,ny;
        std::cout << "Reading nx and ny" << std::endl;
        try{
            nx = std::stoi(first_line[0]);
            ny = std::stoi(first_line[1]);
        } catch(...) {
            return pw::DataSignature::UNKNOWN;
        }
        std::cout << "nx " << std::to_string(nx) << std::endl;
        std::cout << "ny " << std::to_string(ny) << std::endl;
        std::vector<std::string> line;
        for(auto i = 0; i < (nx+ny-2); i++){
            getLineOfData(fin,line);
            if(line.size() != 1)
                return pw::DataSignature::UNKNOWN;
        }
        getLineOfData(fin,line);
        for(const auto& val : line)
            std::cout << "data line value: " << val << std::endl;
        std::cout << "line size: " << line.size() << std::endl;
        if(line.size() < 2)
            return pw::DataSignature::UNKNOWN;
        else if(line.size() == ny)
            return pw::DataSignature::XYZ;
        else if(line.size() == 2*ny)
            return pw::DataSignature::XYCVZ;
        else
            return pw::DataSignature::UNKNOWN;
    }
    return pw::DataSignature::UNKNOWN;
}


pw::OperatorSignature operatorSignature(const std::filesystem::path& path)
{
    std::ifstream stream{path};
    pw::metadataMap meta_map = getHeaderContent(stream);
    if(meta_map.find("OperatorSignature") != meta_map.end())
        return static_cast<pw::OperatorSignature>(std::stoi(meta_map["OperatorSignature"]));
    else
        return pw::OperatorSignature::NONE;
}

pw::metadataMap getHeaderContent(std::ifstream& fin){
    std::map<std::string,std::string> header_map;
    std::string line_feed;
    int oldpos = fin.tellg();
    while(std::getline(fin,line_feed)){
       line_feed = pw::eatWhiteSpace(line_feed);
       if(line_feed.empty()){
           oldpos = fin.tellg();
       } else if(line_feed[0] == '#') {
           int char_count = pw::countCharacters(line_feed,':');
           if(char_count != 1)
               continue;
           std::vector<std::string> parsed_param = pw::parseString(line_feed,':');
           std::string param_name = pw::eatWhiteSpace(parsed_param[0]," \t#");
           std::string param_val = pw::eatWhiteSpace(parsed_param[1]);
           header_map[param_name] = param_val;
           oldpos = fin.tellg();
       } else{
           fin.seekg(oldpos);
           return header_map;
       }
    }
    return header_map;
}

void getLineOfData(std::ifstream& fin,std::vector<std::string>& line_data)
{
    std::string line_feed;
    while(std::getline(fin,line_feed)){
       line_feed = pw::eatWhiteSpace(line_feed);
       // Ignore empty lines and comment lines
       if(line_feed.empty() || line_feed[0] == '#')
           continue;
       else{
           // BINGO, found data line
           line_data = pw::parseString(line_feed,' ');
           return;
       }
    }
    line_data.clear();
}

    

}







