
#ifndef SPIDA_SHAPET_H_
#define SPIDA_SHAPET_H_

#include <string>
#include <complex>
#include "spida/shape/shape.hpp"
#include "spida/grid/gridT.h"

namespace spida{

class ShapeT : public Shape1D<double,std::complex<double>>
{
    public:
        ShapeT(double A,double tp,double omega0) :
            Shape1D<double,std::complex<double>>(),m_A(A),m_tp(tp),m_omega0(omega0),
                m_slow_phase(0.0),m_offset(0.0),m_chirp(0.0) {}
        virtual ~ShapeT() {};
        void setChirp(double v) {m_chirp = v;}
        void setOffset(double v) {m_offset = v;}
        void setSlowPhase(double v) {m_slow_phase = v;}
        double amplitude() const {return m_A;}
        double width() const {return m_tp;}
        double offset() const {return m_offset;}
        double chirp() const {return m_chirp;}
        double fastPhase() const {return m_omega0;}
        double slowPhase() const {return m_slow_phase;}

        virtual std::complex<double> compute(double t) const = 0;
        dcmplx envelope(double t) const {return compute(t)*exp(ii*m_omega0*(t-m_offset));}
        void envelope(const std::vector<double>& t,std::vector<dcmplx>& y) const;
        double computeReal(double t) const {return compute(t).real();}
        void computeReal(const std::vector<double>& t,std::vector<double>& y) const;
        dcmplx computePhaseFactor(double t) const;

        void compute(const GridT& grid,std::vector<dcmplx>& y){
            Shape1D::compute(grid.getT(),y);
        }
        void computeReal(const GridT& grid,std::vector<double>& y){
            computeReal(grid.getT(),y);
        }

    private:
        double m_A;
        double m_tp;
        double m_omega0; 
        double m_slow_phase; 
        double m_offset;
        double m_chirp;
};

class GaussT : public ShapeT
{
    public:
        GaussT(double A,double tp,double omega0) : ShapeT(A,tp,omega0) {}
        ~GaussT() {}; 
        std::complex<double> compute(double t) const;
};


class SechT : public ShapeT
{
    public:
        SechT(double A,double tp,double omega0) : ShapeT(A,tp,omega0) {}
        ~SechT() {}; 
        std::complex<double> compute(double t) const;
};

class AiryT : public ShapeT
{
    public:
        AiryT(double A,double tp,double omega0) : ShapeT(A,tp,omega0) {m_apod = 0.25;}
        ~AiryT() {}; 
        std::complex<double> compute(double t) const;
        void setApodization(double val) {m_apod = val;}
        double apodization() const {return m_apod;}
    private:
        double m_apod;
};

class SuperGaussT : public ShapeT
{
    public:
        SuperGaussT(double A,double tp,double omega0,double m) : ShapeT(A,tp,omega0) {m = m_M;}
        ~SuperGaussT() {}; 
        std::complex<double> compute(double t) const;
        void setM(double val) {m_M = val;}
    private:
        int m_M;

};

class BesselT : public ShapeT
{
    public:
        BesselT(double A,double tp,double omega0);
        ~BesselT() {}; 
        std::complex<double> compute(double t) const;
        void setApodization(double val) {m_apod = val;}
        double apodization() const {return m_apod;}
    private:
        double m_j1;
        double m_apod;
};



}

#endif





