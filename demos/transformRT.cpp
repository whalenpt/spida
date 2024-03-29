/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/17/21
 *    Description: Examples of HankelFFTRRVTRT transform
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/besselR.h>
#include <spida/grid/uniformRVT.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelfftRRVT.h>
#include <pwutils/report/dat.hpp>
#include <pwutils/report/reporthelper.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>

//------------------------------------------------------------------------------


int main()
{

    using spida::dcmplx;
    int nr = 200;
    int nt = 256;
    double w0 = 20.0e-6;
    double I0 = 5.0e16;
    double tp = 5.0e-15;
    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridRVT gridT(nt,minT,maxT,minST,maxST);
    spida::BesselRootGridR gridR(nr,6*w0);

    spida::GaussT shapeT(gridT,std::sqrt(I0),tp);
    shapeT.setFastPhase(omega0);
    spida::GaussR shapeR(gridR,1.0,w0);

    auto u0t = shapeT.shapeRV();
    auto u0r = shapeR.shapeRV();
    std::vector<double> u(nr*nt);

    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt + j] = u0r[i]*u0t[j];

    // Report inputs
    dat::ReportData2D in_report{"RT",gridR.getR(),gridT.getT(),u};
    in_report.setItem("xlabel","r");
    in_report.setItem("ylabel","t");

    std::ofstream os;
    os << std::scientific << std::setprecision(3);
    std::cout << in_report.path() << std::endl;
    os << in_report;

    dat::ReportData1D r_rpd{"R",gridR.getR(),u0r};
    dat::ReportData1D t_rpd{"T",gridT.getT(),u0t};
    os << r_rpd;
    os << t_rpd;


    // Compute transforms
    int nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<double> uop(nr*nt);

    spida::HankelFFTRRVT transform(gridR,gridT);
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,uop);

    dat::ReportComplexData2D<double,double,double> out_report("SRST",gridR.getSR(),gridT.getST(),v);
    out_report.setItem("xlabel","kr");
    out_report.setItem("ylabel","omega");

    os << std::scientific << std::setprecision(8);
    std::cout << out_report.path() << std::endl;
    os << out_report;

    dat::ReportData2D rinop{"RTb",gridR.getR(),gridT.getT(),uop};
    rinop.setItem("xlabel","r");
    rinop.setItem("ylabel","t");

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

    dat::ReportComplexData2D rst_report{"RST",gridR.getR(),gridT.getST(),w};
    rst_report.setItem("xlabel","r");
    rst_report.setItem("ylabel","omega");
    std::cout << rst_report.path() << std::endl;
    os << rst_report;

    dat::ReportData2D srt_report{"SRT",gridR.getSR(),gridT.getT(),zeta};

    srt_report.setItem("xlabel","kr");
    srt_report.setItem("ylabel","t");
    std::cout << srt_report.path() << std::endl;
    os << srt_report;

    dat::ReportComplexData2D rst_report2{"RSTb",gridR.getR(),gridT.getST(),w2};
    rst_report2.setItem("xlabel","r");
    rst_report2.setItem("ylabel","omega");
    os << rst_report2;

    dat::ReportData2D srt_report2{"SRTb",gridR.getSR(),gridT.getT(),zeta2};
    srt_report2.setItem("xlabel","kr");
    srt_report2.setItem("ylabel","t");
    os << srt_report2;

    // Multithread speedup
    unsigned MAX_THREADS = 8;
    unsigned NUM_LOOPS = 20;
    std::vector<int> rt_srst_timings(MAX_THREADS);
    std::vector<int> srst_rt_timings(MAX_THREADS);
    for(unsigned threads = 1; threads < MAX_THREADS; threads++){
        spida::HankelFFTRRVT transform_threaded(gridR,gridT,threads);

        auto start_time = std::chrono::steady_clock::now();
        for(unsigned i = 0; i < NUM_LOOPS; i++)
            transform_threaded.RT_To_SRST(u,v);
        auto elapsed_time = std::chrono::steady_clock::now() - start_time;
        auto rt_srst_timing = std::chrono::duration_cast<std::chrono::microseconds>(elapsed_time).count();
        rt_srst_timings[threads] = static_cast<int>(rt_srst_timing);

        start_time = std::chrono::steady_clock::now();
        for(unsigned i = 0; i < NUM_LOOPS; i++)
            transform_threaded.SRST_To_RT(v,u);
        elapsed_time = std::chrono::steady_clock::now() - start_time;

        auto srst_rt_timing = std::chrono::duration_cast<std::chrono::microseconds>(elapsed_time).count();
        srst_rt_timings[threads] = static_cast<int>(srst_rt_timing);
    }

    for(unsigned threads = 1; threads < MAX_THREADS; threads++){
        std::cout << "HankelFFTRRVT RT_To_SRST duration with " << threads << " thread(s): "\
                  << rt_srst_timings[threads] << "us" << std::endl;
    }
    std::cout << std::endl;
    for(unsigned threads = 1; threads < MAX_THREADS; threads++){
        std::cout << "HankelFFTRRVT SRST_To_RT duration with " << threads << " thread(s): "\
                  << srst_rt_timings[threads] << "us" << std::endl;
    }
    std::cout << std::endl;

    return 0;
}









