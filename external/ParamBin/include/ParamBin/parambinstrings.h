
#ifndef PARAMBINSTRINGS_H_ 
#define PARAMBINSTRINGS_H_

#include<vector>
#include<string>
#include<iostream>

namespace pwbin{

//  using namespace std;
  std::string trimString(const std::string& str,std::string delim = " \t");
  int countFirstChar(const std::string& str,std::string delim = " \t");

  std::string decommentString(const std::string& str,std::string delim = "#");
  std::string eatWhiteSpace(const std::string& str,std::string delim = " \t",
    std::string fill = " ");

  std::vector<std::string> parseString(const std::string& str,const char delim); 
  std::string joinVector(const std::vector<std::string>& vec,const char delim);

  bool isWhitespace(const std::string& str);
  std::vector<std::string> subStrings(const std::string& str,const char delim1 = '(',const char delim2 = ')');
  std::string subString(const std::string& str,const char delim1 = '(',const char delim2 = ')'); 

  bool findString(std::string str,std::string str2); 
  bool splitString(std::string str,std::string str2,std::string& first,std::string& last); 
  bool replaceString(std::string& str,std::string str2,std::string str3); 

  std::string removeSubStrings(const std::string& str,const char delim1 = '(',const char delim2 = ')'); 
  std::string removeSubString(const std::string& str,const char delim1 = '(',const char delim2 = ')');
  std::string parseFirst(std::string& str,const char delim);
  std::string parseLast(std::string& str,const char delim);
  std::string stripLast(const std::string& str,const char delim);
  std::string stripFirst(const std::string& str,const char delim);

  int countWords(const std::string& str);
  int countCharacters(const std::string& str,char c);
  std::string stringLowerCase(const std::string& str);
  std::string stringUpperCase(const std::string& str);

  void addLineBreaks(std::string& str,int charsPerLine);
  std::ostream& columnOutput(std::ostream&,std::string,int,int,int);
}

#endif

