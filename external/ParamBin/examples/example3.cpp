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
#include <ParamBin/parambin.hpp>

//------------------------------------------------------------------------------

int main(int argc,char* argv[])
{

  ParamBin input;
  input.set("var1",1.3);
  input.set("var2",10);
  input.set("var3","string");
  input.setBool("var4",true);
  input.set("var5",1.8);

  char val = 'a';
  std::string val2 = "second value";
  ParamBin group1;
  group1.set("var1",val);
  group1.set("var2",val2);
  input.set("GROUP1",group1);

  ParamBin group2;
  group2.set("var1",5.3e8);
  group2.set("var2","1,2,3");

  ParamBin group3;
  group3.set("var1",1);
  group3.set("var2","subsubgroup");

  group2.set("GROUP3",group3);
  input.set("GROUP2",group2);

  std::vector<std::string> str_vals;
  str_vals.push_back("first val");
  str_vals.push_back("second val");
  str_vals.push_back("third val");
  input.set("string vector",str_vals);

  std::vector<double> dbl_vals;
  dbl_vals.push_back(1.0e-4);
  dbl_vals.push_back(1e4);
  dbl_vals.push_back(3.14);
  input.set("double vector",dbl_vals);
  input.set("int vector","1,2,3");

  std::cout << input << std::endl << std::endl;

  std::cout << "var1 = " << input.getDbl("var1") << std::endl;
  std::cout << "var5 = " << input.getDbl("var5") << std::endl;
  std::cout << "var1 + var5 = " << input.getDbl("var1")+input.getDbl("var5") << std::endl;

  std::cout << "var2 = " << input.getInt("var2") << std::endl;
  std::cout << "var3 = " << input.getStr("var3") << std::endl;
  std::cout << "Uppercase var3 = " << input.getStrU("var3") << std::endl;
  if(input.getBool("var4"))
      std::cout << "var4 is true!" << std::endl;
  else
      std::cout << "var4 is false!" << std::endl;

  std::cout << std::endl;

  std::vector<double> dbl_vec = input.getDblVec("double vector");
  for(int i = 0; i < dbl_vec.size(); i++)
      std::cout << "double vector[" << i << "] = " << dbl_vec[i] << std::endl;

  std::cout << std::endl;

  std::vector<std::string> str_vec = input.getStrVec("string vector");
  for(int i = 0; i < str_vec.size(); i++)
      std::cout << "string vector[" << i << "] = " << str_vec[i] << std::endl;

  std::cout << std::endl;

  std::vector<int> int_vec = input.getIntVec("int vector");
  for(int i = 0; i < int_vec.size(); i++)
      std::cout << "int vector[" << i << "] = " << int_vec[i] << std::endl;

  std::cout << std::endl << std::endl;


  return 0;
}







