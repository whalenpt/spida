
#ifndef SPIDA_SHAPEX_H_
#define SPIDA_SHAPEX_H_

#include <string>
#include <complex>
#include "spida/shape/shape.hpp"
#include "spida/grid/gridX.h"

namespace spida{

// ShapeX in form of A0*shape(c0*(x-offset)+phi0)
class ShapeX : public Shape1D<double,double>
{
    public:
        ShapeX(double A0,double c0) :
            Shape1D<double,double>(),m_A0(A0),m_c0(c0)
            ,m_phi0(0.0), m_offset(0.0)
        virtual ~ShapeX() {};
        void setAmplitude(double A0) {m_A0 = A0;}
        void setWidth(double c0) {m_c0 = c0;}
        void setOffset(double offset) {m_offset = offset;}
        void setPhase(double phi0) {m_phi0 = phi0;}
        double amplitude() const {return m_A;}
        double width() const {return m_c0;}
        double phase() const {return m_phi0;}
        double offset() const {return m_offset;}
        virtual double compute(double x) const = 0;
    private:
        double m_A0;
        double m_c0;
        double m_phi0; 
        double m_offset;
};

class ExpX : public ShapeX
{
    public:
        ExpX(double A0,double c0) : ShapeX(A0,c0) {}
        ~ExpX() {}; 
        double compute(double x) const;
};

class SinX : public ShapeX
{
    public:
        SinX(double A0,double c0) : ShapeX(A0,c0) {}
        ~SinX() {}; 
        double compute(double x) const;
};

class CosX : public ShapeX
{
    public:
        CosX(double A0,double c0) : ShapeX(A0,c0) {}
        ~CosX() {}; 
        double compute(double x) const;
};



}

#endif





