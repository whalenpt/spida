/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/17/21
 *    Description: Examples of PeriodicT transform
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC.

#include <iomanip>
#include <iostream>

#include <spida/grid/uniformRVT.h>
#include <spida/shape/shapeT.h>
#include <spida/transform/fftRVT.h>
#include <utils/report.hpp>

//------------------------------------------------------------------------------


int main()
{

    using spida::dcmplx;
    unsigned nt = 1024;
    double I0 = 5.0e16;
    double tp = 5.0e-15;
    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridRVT gridT(nt,minT,maxT,minST,maxST);
    spida::FFTRVT transform(gridT);
    spida::GaussT shapeT(gridT,std::sqrt(I0),tp);
    shapeT.setFastPhase(omega0);

    auto u = shapeT.shapeRV();
    unsigned nst = gridT.getNst();
    std::vector<dcmplx> v(nst);

    transform.T_To_ST(u,v);

    spida::Report1D in_report{"T", gridT.getT(), u};
    spida::ReportComplex1D out_report{"ST", gridT.getST(), v};

    std::ofstream os;
    std::cout << in_report.path() << std::endl;
    os << in_report;

    std::cout << out_report.path() << std::endl;
    os << out_report;

    return 0;
}



