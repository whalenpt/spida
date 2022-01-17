
#include "pwutils/pwmath.hpp"
#include "pwutils/pwstrings.h"
#include <algorithm>
#include <string>
#include <cmath>
#include <stdexcept>

namespace pw{

int intceil(int x,int y)
{
    int q = x/y + (x % y != 0);
    return q;
}

unsigned int intceil(unsigned int x,unsigned int y)
{
    unsigned int q = x/y + (x % y != 0);
    return q;
}

int factorial(int n) 
{
  return (n == 1 || n ==0 ) ? 1 : factorial(n-1)*n;
}

bool isInteger(const std::string& s)
{
    try{
        [[maybe_unused]] int val = std::stoi(s);
        return true;
    }
    catch(std::exception& e){
        return false;
    }
    return true;
}

bool rowIsIntegers(const std::vector<std::string>& row){
    if(row.empty())
        return false;
    for(const auto& item : row)
        if(!isInteger(item))
            return false;
    return true;
}

bool lineIsIntegers(const std::string& s)
{
    std::string line = pw::eatWhiteSpace(s);
    if(line.empty())
        return false;
    std::vector<std::string> line_data = pw::parseString(line,' ');
    for(const auto& item : line_data)
        if(!isInteger(item))
            return false;
    return true;
}

bool isDouble(const std::string& s){
    try{
        [[maybe_unused]] double val = std::stod(s);
        return true;
    }
    catch(std::exception& e){
        return false;
    }
    return true;
}

bool rowIsDoubles(const std::vector<std::string>& row){
    if(row.empty())
        return false;
    for(const auto& item : row)
        if(!isDouble(item))
            return false;
    return true;
}

bool lineIsDoubles(const std::string& s)
{
    std::string line = pw::eatWhiteSpace(s);
    if(line.empty())
        return false;
    std::vector<std::string> line_data = pw::parseString(line,' ');
    for(const auto& item : line_data)
        if(!isDouble(item))
            return false;
    return true;
}




}



