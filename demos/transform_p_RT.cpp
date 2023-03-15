/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 12/28/21
 *    Description: Demo running multithreaded HankelTransform
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
#include <thread>

//------------------------------------------------------------------------------

const int THREADS = std::max(1,static_cast<int>(std::thread::hardware_concurrency()));
const int NUM_LOOPS = 50;

int main()
{

    using spida::dcmplx;
    int nr = 1600;
    int nt = 16384;
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

    // Compute transforms
    int nst = gridT.getNst();

    std::vector<dcmplx> v(nr*nst);
    std::vector<double> uop(nr*nt);

    spida::HankelFFTRRVT transform(gridR,gridT,THREADS);
    for(int i = 0; i < NUM_LOOPS; i++) { 
        transform.RT_To_SRST(u,v);
        transform.SRST_To_RT(v,uop);
    }

    return 0;
}
