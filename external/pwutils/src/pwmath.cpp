
#include "pwutils/pwmath.hpp"
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>

namespace pw{

int factorial(int n) 
{
  return (n == 1 || n ==0 ) ? 1 : factorial(n-1)*n;
}

bool isInteger(const std::string& s)
{
    return (s.find_first_not_of("0123456789") == std::string::npos);
}

bool rowIsIntegers(const std::vector<std::string>& row){
    if(row.empty())
        return false;
    for(unsigned int i = 0; i < row.size(); i++){
        if(!isInteger(row[i]))
            return false;
    }
    return true;
}

bool isDouble(const std::string& s){
    std::istringstream iss(s);
    double d;
    return iss >> d && !iss.ignore();
}

bool rowIsDoubles(const std::vector<std::string>& row){
    if(row.empty())
        return false;
    for(unsigned int i = 0; i < row.size(); i++){
        if(!isDouble(row[i]))
            return false;
    }
    return true;
}



}


