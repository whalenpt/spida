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
#include <spida/grid/uniformT.h>
#include <iostream>

//------------------------------------------------------------------------------


int main()
{
    std::cout << std::endl << "BesselRootGridR..." << std::endl;
    int N = 8;
    spida::BesselRootGridR grid(N,1.0);
    std::cout << "Bessel Root Grid: " << std::endl;
    for(const auto& item : grid.getR())
        std::cout << item << std::endl;
    std::cout << "Bessel Root Grid Zeros: " << std::endl;
    for(const auto& item : grid.getBesselRoots())
        std::cout << item << std::endl;


    std::cout << std::endl << "UniformGridT..." << std::endl;
    int nt = 1024;
    double tp = 5.0e-15;
//    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridT gridT(nt,minT,maxT,minST,maxST);
    const std::vector<double>& t = gridT.getT();
    const std::vector<double>& st = gridT.getST();
    std::cout << "MinT: " << gridT.getMinT() << std::endl;
    std::cout << "MaxT: " << gridT.getMaxT() << std::endl;
    std::cout << "MinT: " << t[0] << std::endl;
    std::cout << "MaxT: " << t.back() << std::endl;
    std::cout << "MinST: " << gridT.getMinST() << std::endl;
    std::cout << "MaxST: " << gridT.getMaxST() << std::endl;
    std::cout << "MinST: " << st[0] << std::endl;
    std::cout << "MaxST: " << st.back() << std::endl;
    std::cout << std::endl << "T-grid:" << std::endl;
    for(auto i = 0; i < nt; i+=nt/10)
        std::cout << t[i] << std::endl;
    std::cout << std::endl << "ST-grid:" << std::endl;
    for(auto i = 0; i < st.size(); i+=st.size()/10)
        std::cout << st[i] << std::endl;

    return 0;
}








