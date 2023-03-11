//shapeX.h
#pragma once

#include <complex>
#include <string>
#include "spida/grid/gridX.h"
#include "spida/helper/constants.h"
#include "spida/shape/shape.h"

namespace spida{

class ShapeX : public Shape
{
    public:
        ShapeX(const GridX& grid,double A,double w0); 
        ~ShapeX() override = default;
        void setAmplitude(double A) {m_A = A;}
        void setWidth(double w0) {m_w0 = w0;}
        void setOffset(double offset) {m_offset = offset;}
        void setPhase(double phi0) {m_phi0 = phi0;}
        double amplitude() const override {return m_A;}
        double width() const {return m_w0;}
        double phase() const {return m_phi0;}
        double offset() const {return m_offset;}

        std::vector<dcmplx> shapeCV() const;
        std::vector<double> shapeRV() const;
        const std::vector<double>& getX() const {return m_x;} 

    private:
        std::vector<double> m_x;
        double m_A;
        double m_w0;
        double m_offset{0.0};
        double m_phi0{0.0};
        virtual double compute(double x) const = 0;
};

class ExpX : public ShapeX
{
    public:
        using ShapeX::ShapeX;
        ~ExpX() override = default;
    private:
        double compute(double x) const override;
};

class SinX : public ShapeX
{
    public:
        using ShapeX::ShapeX;
        ~SinX() override = default;
    private:
        double compute(double x) const override;
};

class CosX : public ShapeX
{
    public:
        using ShapeX::ShapeX;
        ~CosX() override = default;
    private:
        double compute(double x) const override;
};

}