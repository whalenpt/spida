/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/17/21
 *    Description: Examples of using Grid classes
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/besselR.h>
#include <iostream>

//------------------------------------------------------------------------------


int main()
{
    int N = 8;
    spida::BesselRootGridR grid(N,1.0);
    std::cout << "Bessel Root Grid: " << std::endl;
    for(const auto& item : grid.getR())
        std::cout << item << std::endl;
    std::cout << "Bessel Root Grid Zeros: " << std::endl;
    for(const auto& item : grid.getBesselRoots())
        std::cout << item << std::endl;

    return 0;
}








