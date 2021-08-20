
#ifndef SPIDA_SHAPER_H_
#define SPIDA_SHAPER_H_

#include <string>
#include <complex>
#include "spida/shape/shape.hpp"
#include "spida/grid/gridR.h"

namespace spida{

class ShapeR;
void compute(const GridR& grid,const ShapeR& shape,std::vector<double>& out);

// ShapeR in form of A0*shape(r/w0)
class ShapeR : public Shape1D<double,double>
{
    public:
        ShapeR(double A0,double w0,double offset) :
            Shape1D<double,double>(),m_A0(A0),m_w0(w0),m_offset(offset) {}
        virtual ~ShapeR() {};
        void setAmplitude(double A0) {m_A0 = A0;}
        void setWidth(double w0) {m_w0 = w0;}
        void setOffset(double offset) {m_offset = offset;}
        double amplitude() const {return m_A0;}
        double width() const {return m_w0;}
        double offset() const {return m_offset;}
        virtual double compute(double r) const = 0;
        void compute(const GridR& grid,std::vector<double>& y) const;
    private:
//        GridR m_grid;
        double m_A0;
        double m_w0;
        double m_offset;
};

class GaussR : public ShapeR
{
    public:
        GaussR(double A,double w0,double offset=0.0) : ShapeR(A,w0,offset) {}
        ~GaussR() {}; 
        double compute(double r) const;
};


}

#endif




