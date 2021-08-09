
#ifndef SPIDA_SHAPE_HPP_
#define SPIDA_SHAPE_HPP_

#include <string>
#include <complex>
#include <vector>
#include <cassert>
#include "spida/constants.h"

namespace spida{

template<class T1,class T2>
class Shape1D 
{
    public:
        Shape1D() {};
        virtual ~Shape1D() {};
        virtual T2 compute(T1 x) const = 0;
        void compute(const std::vector<T1>& x,std::vector<T2>& y) const;
};

template<class T1,class T2>
void Shape1D<T1,T2>::compute(const std::vector<T1>& x,std::vector<T2>& y) const
{
    assert(x.size() == y.size());
    for(auto i = 0; i < x.size(); i++)
        y[i] = compute(x[i]);
}


//class Shape2D
//{
//  public:
//    Shape2D() {};
//    virtual ~Shape2D() {};
//    virtual double compute(double x,double y) const = 0;
//    virtual void compute(const double* x,const double* y,double* z,int nx,int ny) const = 0;
//};

}

#endif





