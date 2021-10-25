/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/22/21
 *    Description: Example of HankelR transform
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/besselR.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelR.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/reporthelper.h>
#include <iostream>
#include <fstream>

//------------------------------------------------------------------------------

int main()
{
    // dcmplx is short for std::complex<double>
    int nr = 200;
    double w0 = 20.0e-6;
    double I0 = 5.0e16;

    // Set up a grid of Bessel function roots
    spida::BesselRootGridR gridR(nr,12*w0);
    // GaussR is a Shape class used for gaussian functions with a Bessel root grid
    spida::GaussR shapeR(gridR,std::sqrt(I0),w0);
    // HankelTransformR is a Hankel transform using Bessel function roots as the grid
    spida::HankelTransformR transform(gridR);

    std::vector<double> u = shapeR.shapeRV();
    std::vector<double> v(nr);
    // Transform to spectral space
    transform.R_To_SR(u,v);

    dat::ReportData1D<double,double> in_report("R",gridR.getR(),u);
    in_report.setItem("xlabel","r");
    dat::ReportData1D<double,double> out_report("SR",gridR.getSR(),v);
    out_report.setItem("xlabel","kr");

    std::ofstream os;
    // can set precision for reporting
    os << std::scientific << std::setprecision(3);
    os << in_report;
    os << std::scientific << std::setprecision(8);
    os << out_report;

    std::cout << in_report.path().string() << std::endl;
    std::cout << out_report.path().string() << std::endl;

    using spida::dcmplx;
    std::vector<dcmplx> ucmplx(nr);
    for(auto i = 0; i < nr; i++)
        ucmplx[i] = dcmplx(u[i],0.0);
    std::vector<dcmplx> vcmplx(nr);

    transform.R_To_SR(ucmplx,vcmplx);
    dat::ReportComplexData1D<double,double> cmplx_report("SRcmplx",gridR.getSR(),vcmplx);
    cmplx_report.setItem("xlabel","kr");
    std::cout << cmplx_report.path() << std::endl;
    os << cmplx_report;
    os.close();

    return 0;
}

