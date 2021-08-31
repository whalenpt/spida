
#ifndef SPIDA_SHAPEX_H_
#define SPIDA_SHAPEX_H_

#include <string>
#include <complex>
#include "spida/shape/shape.h"
#include "spida/grid/gridX.h"

namespace spida{

class ShapeX : public Shape
{
    public:
        ShapeX(const GridX& grid,double A0,double c0) :
            Shape(grid),
            m_grid(grid),
            m_A(A0),
            m_c0(c0), 
            m_offset(0.0), 
            m_phi0(0.0)
        virtual ~ShapeX() {};
        void setAmplitude(double A0) {m_A = A0;}
        void setWidth(double c0) {m_c0 = c0;}
        void setOffset(double offset) {m_offset = offset;}
        void setPhase(double phi0) {m_phi0 = phi0;}
        double amplitude() const {return m_A;}
        double width() const {return m_c0;}
        double phase() const {return m_phi0;}
        double offset() const {return m_offset;}
        void shape(std::vector<double>& v) const;

    private:
        GridX m_grid;
        double m_A0;
        double m_c0;
        double m_offset;
        double m_phi0; 
        virtual double compute(double x) const = 0;
};

class ExpX : public ShapeX
{
    public:
        ExpX(double A0,double c0) : ShapeX(A0,c0) {}
        ~ExpX() {}; 
    private:
        double compute(double x) const;
};

class SinX : public ShapeX
{
    public:
        SinX(double A0,double c0) : ShapeX(A0,c0) {}
        ~SinX() {}; 
    private:
        double compute(double x) const;
};

class CosX : public ShapeX
{
    public:
        CosX(double A0,double c0) : ShapeX(A0,c0) {}
        ~CosX() {}; 
    private:
        double compute(double x) const;
};



}

#endif





