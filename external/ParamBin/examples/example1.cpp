/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 04/09/2021
 *    Description: "Reads a ParamBin file."
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <iostream>
#include <vector>
#include <string>
#include <ParamBin/parambin.hpp>

//------------------------------------------------------------------------------

int main(int argc,char* argv[])
{
  if(argc > 1){
    std::cout << "Filename: " << argv[1] << std::endl;
    try{
      ParamBin input(argv[1]);
			std::cout << "Here is your ParamBin Printed!!!" << std::endl << std::endl;
			std::cout << input << std::endl;
		}
		catch(ParamBinException& e){
      std::cout << "A ParamBinException was caught. Here's the message:" << std::endl;
      std::cout << e.what() << std::endl;
		} catch(...) {
			std::cout << "An exception was caught. Failed to read the file." << std::endl;
		}
  }
  else 
  {
    std::cerr << "Please specify an input file " << std::endl;
    std::cerr << "ABORTING" << std::endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}







