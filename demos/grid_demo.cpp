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
#include <spida/grid/uniformRVT.h>
#include <spdlog/spdlog.h>

//------------------------------------------------------------------------------


int main()
{
    spdlog::info("");
    spdlog::info("BesselRootGridR...");
    unsigned N = 8;
    spida::BesselRootGridR grid(N,1.0);
    spdlog::info("Bessel Root Grid: ");
    for(const auto& item : grid.getR())
        spdlog::info("{}", item);
    spdlog::info("Bessel Root Grid Zeros: ");
    for(const auto& item : grid.getBesselRoots())
        spdlog::info("{}", item);


    spdlog::info("");
    spdlog::info("UniformGridT...");
    unsigned nt = 1024;
    double tp = 5.0e-15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridRVT gridT(nt,minT,maxT,minST,maxST);
    const std::vector<double>& t = gridT.getT();
    const std::vector<double>& st = gridT.getST();
    spdlog::info("MinT: {}", gridT.getMinT());
    spdlog::info("MaxT: {}", gridT.getMaxT());
    spdlog::info("MinT: {}", t[0]);
    spdlog::info("MaxT: {}", t.back());
    spdlog::info("MinST: {}", gridT.getMinST());
    spdlog::info("MaxST: {}", gridT.getMaxST());
    spdlog::info("MinST: {}", st[0]);
    spdlog::info("MaxST: {}", st.back());
    spdlog::info("");
    spdlog::info("T-grid:");
    for(unsigned i = 0; i < nt; i+=nt/10)
        spdlog::info("{}", t[i]);
    spdlog::info("");
    spdlog::info("ST-grid:");
    for(unsigned i = 0; i < st.size(); i+=st.size()/10)
        spdlog::info("{}", st[i]);

    return 0;
}