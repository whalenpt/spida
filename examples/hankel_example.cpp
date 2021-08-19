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
#include <spida/transform/hankelR.h>
#include <cmath>
#include <iostream>

//------------------------------------------------------------------------------


int main()
{
    int N = 20;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);
    //spida::printHankel(transform);

    std::vector<double> in(N);
    std::vector<double> in2(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::cout << "jN : " << grid.getjN() << std::endl;
    std::cout << "maxR : " << grid.getMaxR() << std::endl;
    std::cout << "maxSR : " << grid.getMaxSR() << std::endl;

    double a = 5.0;
    for(auto i = 0; i < N; i++)
        in[i] = exp(-pow(a*r[i],2));
    for(auto i = 0; i < N; i++){
        double beta = 1.0/(2.0*pow(a,2));
        exact[i] = beta*exp(-pow(kr[i],2)/(4.0*pow(a,2)));
    }
    transform.R_To_SR(in,out);
    transform.SR_To_R(out,in2);

    for(auto i = 0; i < N; i++)
        std::cout << "OUT: " << out[i] << " - " << "EXPECT: " << exact[i] << std::endl;

    for(auto i = 0; i < N; i++)
        std::cout << "IN: " << in[i] << " - " << "IN2: " << in2[i] << std::endl;

    N = 256;
    rmax = 32.0;
    a = 5.0;

    in.clear();
    out.clear();
    exact.clear();

    in.resize(N);
    out.resize(N);
    exact.resize(N);

    spida::BesselRootGridR grid2(N,rmax);
    spida::HankelTransformR transform2(grid2);

    const std::vector<double>& r2 = grid2.getR();
    const std::vector<double>& kr2 = grid2.getSR();

    for(auto i = 0; i < N; i++)
        in[i] = r2[i] < 1e-4 ? 1.0-pow(a*r2[i],2) : sin(a*r2[i])/(a*r2[i]);
    for(auto i = 0; i < N; i++)
        exact[i] = kr2[i] < a ? 1.0/(pow(a,2)*sqrt(1.0-pow(kr2[i]/a,2))) : 0.0;
    transform2.R_To_SR(in,out);

    for(auto i = 0; i < N; i++)
        std::cout << "OUT: " << out[i] << " - " << "EXPECT: " << exact[i] << std::endl;

    return 0;
}









