
#ifndef SPIDA_SHAPET_H_
#define SPIDA_SHAPET_H_

#include <string>
#include <complex>
#include "spida/helper/constants.h"
#include "spida/grid/gridT.h"

namespace spida{

// Interface class
class ShapeT 
{
    public:
        ShapeT(const GridT& grid,double A,double tp);
        virtual ~ShapeT() {};
        void setChirp(double v) {m_chirp = v;}
        void setOffset(double v) {m_offset = v;}
        void setSlowPhase(double v) {m_slow_phase = v;}
        void setAmplitude(double v) {m_A = v;}
        void setWidth(double v) {m_tp = v;}
        void setFastPhase(double v) {m_omega0 = v;}

        double amplitude() const {return m_A;}
        double width() const {return m_tp;}
        double offset() const {return m_offset;}
        double chirp() const {return m_chirp;}
        double fastPhase() const {return m_omega0;}
        double slowPhase() const {return m_slow_phase;}

        void shape(std::vector<dcmplx>& v) const;
        void shapeReal(std::vector<double>& v) const;
        void envelope(std::vector<dcmplx>& v) const;

    private:
        std::vector<double> m_t;
        double m_A;
        double m_tp;
        double m_offset;
        double m_chirp;
        double m_slow_phase; 
        double m_omega0; 

        dcmplx fastPhaseFactor(double t) const; 
        dcmplx slowPhaseFactor(double t) const; 

        dcmplx computeShape(double t) const {return computeEnvelope(t)*fastPhaseFactor(t);}
        double computeShapeReal(double t) const {return computeShape(t).real();}
        dcmplx computeEnvelope(double t) const {return compute(t)*slowPhaseFactor(t);}

        // compute, will compute base shape
        virtual double compute(double t) const = 0;
};

class GaussT : public ShapeT
{
    public:
        GaussT(const GridT& grid,double A,double tp) : 
            ShapeT(grid,A,tp) {}
        ~GaussT() {}; 
    private:
        double compute(double t) const;
};


class SechT : public ShapeT
{
    public:
        SechT(const GridT& grid,double A,double tp) : 
            ShapeT(grid,A,tp) {}
        ~SechT() {}; 
    private:
        double compute(double t) const;
};

class AiryT : public ShapeT
{
    public:
        AiryT(const GridT& grid,double A,double tp,double apod) : 
            ShapeT(grid,A,tp), m_apod(apod) {}
        ~AiryT() {}; 
        void setApodization(double val) {m_apod = val;}
        double apodization() const {return m_apod;}
    private:
        double compute(double t) const;
        double m_apod;
};

class SuperGaussT : public ShapeT
{
    public:
        SuperGaussT(const GridT& grid,double A,double tp,double m) : 
            ShapeT(grid,A,tp), m_M(m) {}
        ~SuperGaussT() {}; 
        void setM(double val) {m_M = val;}
    private:
        double compute(double t) const;
        int m_M;

};

class BesselT : public ShapeT
{
    public:
        BesselT(const GridT& grid,double A,double tp,double apod);
        ~BesselT() {}; 
        void setApodization(double val) {m_apod = val;}
        double apodization() const {return m_apod;}
    private:
        double compute(double t) const;
        double m_apod;
        double m_j1;
};



}

#endif





