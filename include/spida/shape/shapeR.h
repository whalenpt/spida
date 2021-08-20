
#ifndef SPIDA_SHAPER_H_
#define SPIDA_SHAPER_H_

#include <string>
#include <complex>
#include "spida/grid/gridR.h"

namespace spida{

// ShapeR in form of A0*shape(r/w0)
class ShapeR
{
    public:
        ShapeR(const GridR& grid,double A,double w0); 
        virtual ~ShapeR() {};
        void setAmplitude(double A) {m_A = A;}
        void setWidth(double w0) {m_w0 = w0;}
        double amplitude() const {return m_A;}
        double width() const {return m_w0;}
        void shape(std::vector<double>& v) const;

    private:
        std::vector<double> m_r;
        double m_A;
        double m_w0;
        virtual double compute(double r) const = 0;
};

class GaussR : public ShapeR
{
    public:
        GaussR(const GridR& grid,double A,double w0) : ShapeR(grid,A,w0) {}
        ~GaussR() {}; 
        double compute(double r) const;
};


}

#endif




