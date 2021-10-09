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

#include <spida/grid/uniformRVT.h>
#include <spida/shape/shapeT.h>
#include <spida/transform/fftRVT.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/reporthelper.h>
#include <iostream>

//------------------------------------------------------------------------------


int main()
{

    using spida::dcmplx;
    int nt = 1024;
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

    std::vector<double> u(nt);
    shapeT.shapeRV(u);
    int nst = gridT.getNst();
    std::vector<dcmplx> v(nst);

    transform.T_To_ST(u,v);

    dat::ReportData1D<double,double> in_report("T",gridT.getT(),u);
    dat::ReportComplexData1D<double,double> out_report("ST",gridT.getST(),v);

    std::ofstream os;
    os << std::scientific << std::setprecision(5);
    std::cout << in_report.path() << std::endl;
    os << in_report;

    os << std::scientific << std::setprecision(8);
    std::cout << out_report.path() << std::endl;
    os << out_report;

    return 0;
}



