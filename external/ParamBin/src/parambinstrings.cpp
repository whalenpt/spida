
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include "ParamBin/parambinstrings.h"

namespace pwbin{

// Check if string is whitespace
bool isWhitespace(const std::string& str){
    if(str.find_first_not_of(' ') != std::string::npos)
        return false;
    else
        return true;
}

void replaceLast(std::string& str,std::string delim1,std::string delim2)
{
  size_t pos1 = str.find_last_of(delim1);
  if(pos1 != std::string::npos) 
    str.replace(pos1,1,delim2);
}

void removeChars(std::string& str,std::string delim,std::string fill)
{
  size_t beginSpace = str.find_first_of(delim);
  while(beginSpace != std::string::npos){
    const size_t endSpace = str.find_first_not_of(delim,beginSpace);
    const size_t range = endSpace - beginSpace;
    str.replace(beginSpace,range,fill);
    const size_t newStart = beginSpace + fill.length();
    beginSpace = str.find_first_of(delim,newStart);
  }
  str = eatWhiteSpace(str);
}

std::string stringLowerCase(const std::string& cstr)
{
    std::string str(cstr);
    std::transform(str.begin(),str.end(),str.begin(), ::tolower);
    return str;
}

std::string stringUpperCase(const std::string& cstr)
{
    std::string str(cstr);
    std::transform(str.begin(),str.end(),str.begin(), ::toupper);
    return str;
}

std::string trimString(const std::string& str,std::string delim)
{
  const size_t pos1 = str.find_first_not_of(delim);
  if(pos1 == std::string::npos) return "";
  const size_t pos2 = str.find_last_not_of(delim);
  const size_t range = pos2 - pos1 + 1;
  return str.substr(pos1,range);
}

int countFirstChar(const std::string& str,std::string delim)
{
  const size_t pos1 = str.find_first_not_of(delim);
  if(pos1 == std::string::npos) return 0;
  std::string fststr = str.substr(0,pos1);
  int cnt = fststr.size();
  return cnt;
}

std::string decommentString(const std::string& str,std::string delim)
{
  size_t pos1 = str.find(delim);
  if(pos1 == std::string::npos) return str;
  return str.substr(0,pos1);
}

std::string eatWhiteSpace(const std::string& str,std::string delim,
    std::string fill)
{
  std::string result(trimString(str,delim));
  size_t beginSpace = result.find_first_of(delim);
  while(beginSpace != std::string::npos){
    const size_t endSpace = result.find_first_not_of(delim,beginSpace);
    const size_t range = endSpace - beginSpace;
    result.replace(beginSpace,range,fill);
    const size_t newStart = beginSpace + fill.length();
    beginSpace = result.find_first_of(delim,newStart);
  }
  return result;
}

std::string joinVector(const std::vector<std::string>& vec,const char delim)
{
    if(vec.empty())
        return "";
    else if(vec.size() == 1)
        return vec[0];

    std::string str = vec[0];
    for(int i = 1; i < vec.size(); i++)
        str += delim + vec[i]; 
    return str;
}

std::vector<std::string> parseString(const std::string& mystr,const char delim) 
{
  std::string str(mystr);
  std::vector<std::string> parsedStr;
  if(str.empty()){
    return parsedStr;
  }
  str = eatWhiteSpace(str);
  size_t pos1 = 0; 
  size_t pos2 = str.find_first_of(delim);
  // If the delimiter character is not found, then add the input string to the output vector and return
  if(pos2 == std::string::npos) {
    parsedStr.push_back(str);
    return parsedStr;
  }
  // Delimiter character is found! 
  // Save the substring of the input string from the start to the first delimiter character
  std::string tempStr = str.substr(0,pos2);
  tempStr = eatWhiteSpace(tempStr);
  // Add this substring to the output vector
  parsedStr.push_back(tempStr);
  // Continue searching for delimiter characters
  while(pos2 != std::string::npos){
    pos1 = pos2+1;
    pos2 = str.find_first_of(delim,pos1);
    size_t range = pos2 - pos1;
    tempStr = str.substr(pos1,range);
    tempStr = eatWhiteSpace(tempStr);
    parsedStr.push_back(tempStr);
  }
  return parsedStr;
}

std::string parseFirst(std::string& str,const char delim)
{
  size_t pos2 = str.find_first_of(delim);
  std::string strippedStr("");
  if(pos2 == std::string::npos) {
    strippedStr = str;
    str.clear();
    return strippedStr;
  }
  strippedStr = str.substr(0,pos2);
  size_t range = std::string::npos - pos2 - 1;
  str = str.substr(pos2+1,range);
  return strippedStr;
}

std::string parseLast(std::string& str,const char delim)
{
  size_t pos1 = str.find_last_of(delim);
  std::string strippedStr(str);
  if(pos1 == std::string::npos) {
    str = "";
    return strippedStr;
  }
  size_t range = std::string::npos - pos1 - 1;
  strippedStr = str.substr(pos1+1,range);
  str = str.substr(0,pos1);
  return strippedStr;
}

std::string stripLast(const std::string& str,const char delim)
{
  size_t pos1 = str.find_last_of(delim);
  std::string strippedStr(str);
  if(pos1 == std::string::npos) {
    return strippedStr;
  }
  size_t range = std::string::npos - pos1 - 1;
  strippedStr = str.substr(pos1+1,range);
  return strippedStr;
}

std::string stripFirst(const std::string& str,const char delim)
{
  size_t pos1 = str.find_last_of(delim);
  if(pos1 == std::string::npos) {
    return "";
  }
  size_t range = std::string::npos - pos1 - 1;
  return str.substr(0,pos1);
}




//std::string parseLast(std::string& str,const char delim)
//{
//  std::vector<std::string> parsed(parseString(str,delim)); 
//  str = parsed.back();
//}

std::vector<std::string> subStrings(const std::string& instr,const char delim1,const char delim2) 
{
  std::vector<std::string> subStr(0);
  if(instr.empty())
    return subStr;

  std::string str = eatWhiteSpace(instr);
  size_t pos1 = str.find_first_of(delim1); 
  if(pos1 == std::string::npos)
    return subStr;

  size_t pos2 = str.find_first_of(delim2,pos1+1);
  if(pos2 == std::string::npos || ((pos2+1) == std::string::npos)) 
    return subStr;
  
  size_t range = pos2-pos1-1;
  std::string tempStr = str.substr(pos1+1,range);
  tempStr = eatWhiteSpace(tempStr);
  subStr.push_back(tempStr);
  while(pos2 != std::string::npos && pos1 != std::string::npos){
    pos1 = str.find_first_of(delim1,pos2+1);
    if(pos1 == std::string::npos)
      return subStr;
    pos2 = str.find_first_of(delim2,pos1+1);
    if(pos2 == std::string::npos) 
      return subStr;
    range = pos2 - pos1 - 1;
    tempStr = str.substr(pos1+1,range);
    tempStr = eatWhiteSpace(tempStr);
    subStr.push_back(tempStr);
  }
  return subStr;
}

std::string subString(const std::string& instr,const char delim1,const char delim2) 
{
  std::string subStr;
  if(instr.empty())
    return subStr;

  std::string str = eatWhiteSpace(instr);
  size_t pos1 = str.find_first_of(delim1); 
  size_t pos2 = str.find_first_of(delim2);
  if(pos1 == std::string::npos || pos2 == std::string::npos) 
    return subStr;
  else{
    size_t range = pos2-pos1-1;
    subStr = str.substr(pos1+1,range);
    subStr = eatWhiteSpace(subStr);
    return subStr;
  }
}

bool findString(std::string str,std::string str2) 
{
  const char delim1 = str2[0];
  const char delim2 = str2[str2.size()-1];
  std::string subStr;
  if(str.empty())
    return false;

  str = eatWhiteSpace(str);
  size_t pos1 = str.find_first_of(delim1); 
  size_t pos2 = str.find_first_of(delim2);
  if(pos1 == std::string::npos || pos2 == std::string::npos) 
    return false;
  size_t range = pos2-pos1+1;
  subStr = str.substr(pos1,range);
  subStr = eatWhiteSpace(subStr);
  if(subStr == str2)
    return true;
  else
    return false;
}

bool splitString(std::string str,std::string str2,std::string& first,std::string& last) 
{
  if(!findString(str,str2)){
    first = str;
    last = "";
    return false;
  }
  const char delim1 = str2[0];
  const char delim2 = str2[str2.size()-1];
  std::string subStr;
  str = eatWhiteSpace(str);
  size_t pos1 = str.find_first_of(delim1); 
  size_t pos2 = str.find_first_of(delim2);
//  size_t range = pos2-pos1+1;
  first = str.substr(0,pos1);
//  std::cout << "First: " << first << std::endl;
  last = str.substr(pos2+1,std::string::npos);
//  std::cout << "Last: " << last << std::endl;
  return true;
}

bool replaceString(std::string& str,std::string str2,std::string str3) 
{
  if(!findString(str,str2))
    return false;

  const char delim1 = str2[0];
  const char delim2 = str2[str2.size()-1];
  std::string subStr;
  str = eatWhiteSpace(str);
  size_t pos1 = str.find_first_of(delim1); 
  size_t pos2 = str.find_first_of(delim2);
  size_t range = pos2-pos1+1;
  subStr = str.substr(pos1,range);
  subStr = eatWhiteSpace(subStr);
  if(subStr == str2){
    str.replace(pos1,range,str3);
    return true;
  }
  else
    return false;
}


std::string removeSubStrings(const std::string& cstr,const char delim1,const char delim2) 
{
  if(cstr.empty())
    return cstr;
  size_t pos1 = cstr.find_first_of(delim1); 
  if(pos1 == std::string::npos)
    return cstr;

  size_t pos2 = cstr.find_first_of(delim2,pos1+1);
  if(pos2 == std::string::npos || ((pos2+1) == std::string::npos))
    return cstr;

  std::string str(cstr);
  while(pos2 != std::string::npos && pos1 != std::string::npos){
    size_t range = pos2-pos1+1;
    str.erase(pos1,range);
    pos1 = str.find_first_of(delim1);
    if(pos1 == std::string::npos)
      return str;
    pos2 = str.find_first_of(delim2,pos1+1);
    if(pos2 == std::string::npos || ((pos2+1) == std::string::npos))
      return str;
  }
  return str;
}

std::string removeSubString(const std::string& cstr,const char delim1,const char delim2) 
{
  if(cstr.empty())
    return cstr;
  std::string str(cstr);
  size_t pos1 = str.find_first_of(delim1); 
  size_t pos2 = str.find_first_of(delim2);
  if(pos2 != std::string::npos && pos1 != std::string::npos){
    size_t range = pos2-pos1+1;
    str.erase(pos1,range);
  }
  return str;
}

int countCharacters(const std::string& str,char c)
{
  int cnt = std::count(str.begin(),str.end(),c);
  return cnt;
}

int countWords(const std::string& str)
{
  int wordCount(0);
  std::stringstream ss(str);
  std::string word;
  while(ss >> word) ++wordCount;
  return wordCount;
}

void addLineBreaks(std::string& str,int charsPerLine)
{
  trimString(str);
  str = eatWhiteSpace(str," \t\n");
  int charNum = charsPerLine;
  int strLength = str.length();
  if(strLength < charsPerLine)
    return;
  while(charNum < strLength){
    str.insert(charNum,1,'\n');
    charNum += charsPerLine;
    strLength++;
  }
}

std::ostream& columnOutput(std::ostream& os,std::string str,int st_col_num,int end_col_num,int init_col_num)
{
  unsigned int width = end_col_num - init_col_num;
  if(str.length() < width){
    os << std::setw(width) << str << std::endl;
    return os;
  }
  std::stringstream ss;
  ss << str;
  unsigned int rowLen;
  std::string tempStr;

  ss >> tempStr;
  tempStr += " ";
  rowLen = tempStr.length();
  os << tempStr;
  tempStr.clear();
  while(rowLen < width){
    ss >> tempStr;
    tempStr += " ";
    rowLen += tempStr.length();
    os << tempStr;
    tempStr.clear();
  }
  os << std::endl;

  width = end_col_num - st_col_num;
  std::string colSpace(st_col_num,' ');
  ss >> tempStr;
  if(ss)
    os << colSpace;
  else 
    return os;
   
  rowLen = 0;
  while(ss) {
    tempStr += " ";
    os << tempStr;
    rowLen += tempStr.length();
    tempStr.clear();
    if(rowLen > width){
      os << std::endl << colSpace;
      rowLen = 0;
    }
    ss >> tempStr;
  }
  os << std::endl;

  return os;
}

std::string NumberToString(int val) {
  std::ostringstream stm;
  stm.precision(12);
  stm << val;
  return stm.str();
}



}


