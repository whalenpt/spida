#pragma once

#include <complex>
#include <stdexcept>
#include <string>
#include "spida/grid/gridR.h"
#include "spida/helper/constants.h"
#include "spida/shape/shape.h"

namespace spida{

// ShapeR in form of A0*shape(r/w0)*exp(-ir^2/f)
//               ->  A0*shape(...)*exp(focusPhaseFactor)
class ShapeR : public Shape
{
    public:
        ShapeR(const GridR& grid,double A,double w0); 
        ~ShapeR() override = default;
        void setAmplitude(double A) {m_A = A;}
        // Precondition: w0 != 0.0  (throws std::domain_error otherwise)
        void setWidth(double w0);
        void setFocus(double f);
        double amplitude() const override {return m_A;}
        double width() const {return m_w0;}
        double focus() const {return m_f;}
        [[nodiscard]] std::vector<dcmplx> shapeCV() const;
        // Returns the real part of shapeCV(). When a focus is set the
        // quadratic phase is non-trivial and the imaginary part is silently
        // discarded; use shapeCV() directly if the full complex field is needed.
        [[nodiscard]] std::vector<double> shapeRV() const;
        const std::vector<double>& getR() const {return m_r;} 

    private:
        std::vector<double> m_r;
        double m_A;
        double m_w0;
        bool m_bool_focus{false};
        double m_f{0.0};
        virtual double compute(double r) const = 0;
        dcmplx focusPhaseFactor(double r) const;
};

class GaussR : public ShapeR
{
    public:
        using ShapeR::ShapeR;
        ~GaussR() override = default;
    private:
        double compute(double r) const override;
};

}
