
#include <cmath>
#include <string>
#include <stdexcept>
#include "spida/grid/chebX.h"
#include "spida/helper/constants.h"

namespace spida{

void setChebExtremaX(double a,double b,std::vector<double>& x)
{
   if(a >= b){
        std::string msg = "Error in setChebExtremaX: interval [a,b] -> " + std::to_string(a) + \
                           " must be less than " + std::to_string(b);
        throw std::invalid_argument(msg);
    }
//    int N = x.size() - 1;
    auto nx = x.size();
    double dc = spida::PI/(nx-1);
    double L = b - a;
    //for(int i = 0; i <= N; i++) x[i] = a + L*(1.0-cos(i*dc))/2.0; 
    for(int j = 0; j < nx; j++) x[nx-j-1] = (cos(j*dc)*L+b+a)/2.0; 
}

void setChebExtremaSX(std::vector<double>& sx)
{
    int nsx = sx.size();
    double dc = spida::PI/(nsx-1);
    for(int j = 0; j < nsx; j++) sx[j] = j*dc;
}

void setChebRootX(double a,double b,std::vector<double>& x)
{
//   if(a >= b){
//        std::string msg = "Error in setChebRootX: interval [a,b] -> a must be less than b.";
//        throw std::invalid_argument(msg);
//    }
   if(a >= b){
        std::string msg = "Error in setChebRootX: interval [a,b] -> " + std::to_string(a) + \
                           " must be less than " + std::to_string(b);
        throw std::invalid_argument(msg);
    }

    auto nx = x.size();
    auto dc = spida::PI/static_cast<double>(nx);
    auto L = b - a;
    //for(int i = 0; i < N; i++) x[i] = a + L*(1.0-cos((i+0.5)*dc))/2.0; 
    for(int j = 0; j < nx; j++) x[nx-j-1] = (cos((j+0.5)*dc)*L+b+a)/2.0; 
}

void setChebRootSX(std::vector<double>& sx)
{
    int nsx = sx.size();
    double dc = spida::PI/static_cast<double>(nsx);
    for(int j = 0; j < nsx; j++) sx[j] = (j+0.5)*dc;
}

ChebExtremaGridX::ChebExtremaGridX(int nx,double min,double max) : ChebGridX(nx,min,max),
    m_x(nx), m_sx(nx)
{
    setChebExtremaX(min,max,m_x);
    setChebExtremaSX(m_sx);
}


ChebRootGridX::ChebRootGridX(int nx,double min,double max) : ChebGridX(nx,min,max),
    m_x(nx), m_sx(nx)
{
    // Grid only includes interior points (not min or max)
    setChebRootX(min,max,m_x);
    setChebRootSX(m_sx);
}

}

