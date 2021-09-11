
#ifndef SPIDA_SHAPER_H_
#define SPIDA_SHAPER_H_

#include <string>
#include <complex>
#include "spida/grid/gridR.h"
#include "spida/shape/shape.h"
#include "spida/helper/constants.h"

namespace spida{

// ShapeR in form of A0*shape(r/w0)*exp(-ir^2/f)
//               ->  A0*shape(...)*exp(focusPhaseFactor)
class ShapeR : public Shape
{
    public:
        ShapeR(const GridR& grid,double A,double w0); 
        virtual ~ShapeR() {};
        void setAmplitude(double A) {m_A = A;}
        void setWidth(double w0) {m_w0 = w0;}
        void setFocus(double f);
        double amplitude() const {return m_A;}
        double width() const {return m_w0;}
        double focus() const {return m_f;}
        std::vector<dcmplx> shapeCV() const;
        std::vector<double> shapeRV() const;
        void shapeCV(std::vector<dcmplx>& v) const;
        void shapeRV(std::vector<double>& v) const;
        const std::vector<double>& getR() const {return m_r;} 

    private:
        std::vector<double> m_r;
        double m_A;
        double m_w0;
        bool m_bool_focus;
        double m_f;
        virtual double compute(double r) const = 0;
        dcmplx focusPhaseFactor(double r) const;
        dcmplx computeShape(double r) const; 
        double computeShapeReal(double r) const {return computeShape(r).real();}
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




