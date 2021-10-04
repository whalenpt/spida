/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/17/21
 *    Description: Examples of using interpolation
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <vector>
#include <cmath>
#include <string>
#include <ParamBin/parambin.hpp>
#include <pwutils/report/dataio.hpp>
#include <spida/helper/interp.h>
//------------------------------------------------------------------------------

int main()
{
    std::vector<double> x(11);
    for(int i = 0; i <=10; i++){
        x[i] = static_cast<double>(i);
    }
    std::vector<double> y(11);
    for(int i = 0; i<= 10; i++){
        y[i] = sin(x[i]);
    }
    spida::LinearInterp interp(x,y);
    int N = 40;
    std::vector<double> xinterp(N);
    for(int i = 0; i < N; i++)
        xinterp[i] = 10.0*i/(N-1);
    std::vector<double> yinterp = interp.eval(xinterp);

    pw::DataIO dataio("outfolder");
    dataio.writeFile("data.dat",x,y);
    dataio.writeFile("interp_data.dat",xinterp,yinterp);

    spida::SplineInterp sinterp(x,y);
    std::vector<double> ysinterp = sinterp.eval(xinterp);
    dataio.writeFile("spline_interp_data.dat",xinterp,ysinterp);

    //Test tridisolve
    int n = 20;
    std::vector<double> a(n-1);
    std::vector<double> b(n);
    std::vector<double> c(n-1);
    std::vector<double> d(n);
    for(int i = 0; i < n-1; i++){
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;
        d[i] = i+1;
    }
    b[n-1] = 2.0;
    d[n-1] = n;
    std::vector<double> r;
    spida::tridisolve(a,b,c,d,r);
    dataio.writeFile("tridisolve.dat",r);

    return 0;
}















