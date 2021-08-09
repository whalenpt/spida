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
  input.set("var5",1.8600);
  input.set("var6",0.000345);
  input.set("var7",100000.00);
  input.set("var8",100005.00001);

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

  std::cout << input << std::endl << std::endl;
  std::cout << "Clearing var1,var3, and var5..." << std::endl;
  input.clear("var1");
  input.clear("var3");
  input.clear("var5");
  std::cout << input << std::endl << std::endl;

  std::cout << "Clearing GROUP1..." << std::endl;
  input.clear("GROUP1");
  std::cout << input << std::endl << std::endl;

  std::cout << "Clearing GROUP3..." << std::endl;
  input.getBin("GROUP2").clear("GROUP3");
  std::cout << input << std::endl << std::endl;

  ParamBin& bin = input.getBin("GROUP2");
  bin.set("new value",2.4445);
  group2.set("other new value","WHERE AM I?");

  std::cout << "ADDED A NEW VALUE! " << std::endl;
  std::cout << input << std::endl << std::endl;

  const ParamBin& const_bin = input.getBin("GROUP2");

  std::cout << "ACCESS A VALUE FROM reference to a constant bin " << std::endl;
  std::cout << "GROUP2 var2 value = " << const_bin.getStr("var2") << std::endl << std::endl;
  std::cout << input << std::endl << std::endl;

//  const_bin.set("var3","value");   //This will fail to compile!


  return 0;
}







