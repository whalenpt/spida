/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/17/21
 *    Description: Examples of HankelPeriodicRT transform
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/besselR.h>
#include <spida/grid/uniformT.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelperiodicRT.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/reporthelper.h>
#include <cmath>
#include <iostream>
#include <fstream>

//------------------------------------------------------------------------------


int main()
{

    using spida::dcmplx;
    int nr = 200;
    int nt = 1024;
    double w0 = 20.0e-6;
    double I0 = 5.0e16;
    double tp = 5.0e-15;
    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridT gridT(nt,minT,maxT,minST,maxST);
    spida::BesselRootGridR gridR(nr,12*w0);

    spida::GaussT shapeT(gridT,std::sqrt(I0),tp);
    shapeT.setFastPhase(omega0);
    spida::GaussR shapeR(gridR,1.0,w0);

    std::vector<double> u0t(nt);
    shapeT.shapeReal(u0t);
    std::vector<double> u0r(nr);
    shapeR.shape(u0r);
    std::vector<double> u(nr*nt);

    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt + j] = u0r[i]*u0t[j];

    // Report inputs
    dat::ReportData2D<double,double,double> in_report("RT",gridR.getR(),gridT.getT(),u);
    in_report.setLabelX("r");
    in_report.setLabelY("t");

    std::ofstream os;
    os << std::scientific << std::setprecision(3);
    std::cout << in_report.path() << std::endl;
    os << in_report;

    dat::ReportData1D<double,double> r_rpd("R",gridR.getR(),u0r);
    dat::ReportData1D<double,double> t_rpd("T",gridT.getT(),u0t);
    os << r_rpd;
    os << t_rpd;


    // Compute transforms
    int nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<double> uop(nr*nt);

    spida::HankelPeriodicTransformRT transform(gridR,gridT);
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,uop);

    dat::ReportComplexData2D<double,double> out_report("SRST",gridR.getSR(),gridT.getST(),v);
    out_report.setLabelX("kr");
    out_report.setLabelY("omega");

    os << std::scientific << std::setprecision(8);
    std::cout << out_report.path() << std::endl;
    os << out_report;

    dat::ReportData2D<double,double,double> rinop("RTb",gridR.getR(),gridT.getT(),uop);
    rinop.setLabelX("r");
    rinop.setLabelY("t");

    os << std::scientific << std::setprecision(3);
    std::cout << rinop.path() << std::endl;
    os << rinop;

    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> w2(nr*nst);
    std::vector<double> zeta(nr*nt);
    std::vector<double> zeta2(nr*nt);

    transform.RT_To_RST(u,w);
    transform.RT_To_SRT(u,zeta);
    transform.SRST_To_RST(v,w2);
    transform.SRST_To_SRT(v,zeta2);


    dat::ReportComplexData2D<double,double> rst_report("RST",gridR.getR(),gridT.getST(),w);
    rst_report.setLabelX("r");
    rst_report.setLabelY("omega");
    std::cout << rst_report.path() << std::endl;
    os << rst_report;

    dat::ReportData2D<double,double,double> srt_report("SRT",gridR.getSR(),gridT.getT(),zeta);
    srt_report.setLabelX("kr");
    srt_report.setLabelY("t");
    std::cout << srt_report.path() << std::endl;
    os << srt_report;

    dat::ReportComplexData2D<double,double> rst_report2("RSTb",gridR.getR(),gridT.getST(),w2);
    rst_report2.setLabelX("r");
    rst_report2.setLabelY("omega");
    os << rst_report2;

    dat::ReportData2D<double,double,double> srt_report2("SRTb",gridR.getSR(),gridT.getT(),zeta2);
    srt_report2.setLabelX("kr");
    srt_report2.setLabelY("t");
    os << srt_report2;

    return 0;
}









