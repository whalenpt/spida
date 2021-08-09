/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 04/09/2021
 *    Description: "Testing of ParamBin class."
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <iostream>
#include <string>
#include <memory>
#include <ParamBin/parambin.hpp>

//------------------------------------------------------------------------------

int main(int argc,char* argv[])
{

  ParamBin bin;
  std::cout << "Alias names can be put in parenthesis." << std::endl;
  // Alias names can be placed in parenthesis
  bin << NamedParam<int>("LongVeryLongLongName (short name)",1000);
  bin << NamedParam<int>("A",10) << NamedParam<int>("B",20) \
      << NamedParam<int>("C",30) << NamedParam<int>("D",40); 
  std::cout << bin << std::endl << std::endl;
  std::cout << "Access LongVeryLongLongName via short name. " << std::endl;
  std::cout << "LongVeryLongLongName = " << bin.getInt("short name") << std::endl;

  ParamBin second_bin;
  second_bin << NamedParam<double>("Var1 (v1)",1.6) << NamedParam<double>("Var2",2.4) \
             << NamedParam<double>("Var3",3.2) << NamedParam<double>("Var4",3.6);

  std::unique_ptr<ParamBin> third_bin(new ParamBin);
  *third_bin << NamedParam<char>("3_1",'X') << NamedParam<char>("3_2",'Y') \
            << NamedParam<char>("3_3 (v3)",'Z'); 


  second_bin.setBin("THIRD BIN",std::move(third_bin));
  bin.setBin("SECOND BIN",second_bin);

  std::cout << bin << std::endl << std::endl;
  // Alias' are searched for recursively in children of the parent,
  // effectively this means the user can gain access to these values from the top of the tree.

  std::cout << "Var1 (v1) = " << bin.getDbl("v1") << std::endl; 
  std::cout << "3_3 (v3) = " << bin.getChar("v3") << std::endl; 


  return 0;
}







