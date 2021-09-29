//shapeX.h
#pragma once

#include <string>
#include <complex>
#include "spida/helper/constants.h"
#include "spida/grid/gridX.h"
#include "spida/shape/shape.h"

namespace spida{

class ShapeX : public Shape
{
    public:
        ShapeX(const GridX& grid,double A,double w0); 
        virtual ~ShapeX() {};
        void setAmplitude(double A) {m_A = A;}
        void setWidth(double w0) {m_w0 = w0;}
        void setOffset(double offset) {m_offset = offset;}
        void setPhase(double phi0) {m_phi0 = phi0;}
        double amplitude() const {return m_A;}
        double width() const {return m_w0;}
        double phase() const {return m_phi0;}
        double offset() const {return m_offset;}

        std::vector<dcmplx> shapeCV() const;
        std::vector<double> shapeRV() const;
        void shapeCV(std::vector<dcmplx>& v) const;
        void shapeRV(std::vector<double>& v) const;
        const std::vector<double>& getX() const {return m_x;} 
        dcmplx shapeCV(double x) const; 
        double shapeRV(double x) const {return shapeCV(x).real();}

    private:
        std::vector<double> m_x;
        double m_A;
        double m_w0;
        double m_offset;
        double m_phi0;
        virtual double compute(double r) const = 0;
};

class ExpX : public ShapeX
{
    public:
        ExpX(const GridX& grid,double A0,double c0) : ShapeX(grid,A0,c0) {}
        ~ExpX() {}; 
    private:
        double compute(double x) const;
};

class SinX : public ShapeX
{
    public:
        SinX(const GridX& grid,double A0,double c0) : ShapeX(grid,A0,c0) {}
        ~SinX() {}; 
    private:
        double compute(double x) const;
};

class CosX : public ShapeX
{
    public:
        CosX(const GridX& grid,double A0,double c0) : ShapeX(grid,A0,c0) {}
        ~CosX() {}; 
    private:
        double compute(double x) const;
};



}






