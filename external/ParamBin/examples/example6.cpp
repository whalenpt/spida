/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
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
  input.set("var1",10);
  int val = input.getInt("var1");
  std::cout << "Int 10: " <<  val << std::endl;

  int intval = 10;
  input.set("var2",intval);
  val = input.getInt("var2");
  std::cout << "Int 10: " <<  val << std::endl;

  return 0;
}







