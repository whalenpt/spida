//ETDCS.h
#pragma once

#include <vector>
#include "spida/rkstiff/solver.h"
#include "spida/helper/constants.h"

namespace spida{

class SolverCS_ETD : public SolverCV_CS
{
    public:
        using SolverCV_CS::SolverCV_CS;
        ~SolverCS_ETD() override = default;
        void setModeCutoff(double val) {m_mode_cutoff = val;}
        void setContourRadius(double val) {m_contour_radi = val;}
        void setContourPoints(int val) {m_contourM = val;}
        double modeCutoff() const {return m_mode_cutoff;}
        double contourRadius() const {return m_contour_radi;}
        int contourPoints() const {return m_contourM;}
    private:
        double m_mode_cutoff{0.01};
        int m_contourM{32};
        double m_contour_radi{1.0};
};

class ETD4 : public SolverCS_ETD
{
    public:
        ETD4(const LinOp& L,const NLfunc& NL,bool use_refs=false);
        ~ETD4() override = default;
    private:
        void updateCoefficients(double dt) noexcept override;
        void updateStages(std::vector<dcmplx>& in) noexcept override;
        int m_sz;
        const std::vector<dcmplx>& L;
        std::vector<dcmplx> EL; std::vector<dcmplx> EL2; 
        std::vector<dcmplx> N1; std::vector<dcmplx> N2; 
        std::vector<dcmplx> N3; std::vector<dcmplx> N4;
        std::vector<dcmplx> tempK;
        std::vector<dcmplx> a21; std::vector<dcmplx> a31; std::vector<dcmplx> a32;
        std::vector<dcmplx> a41; std::vector<dcmplx> a43;
        std::vector<dcmplx> a51; std::vector<dcmplx> a52; std::vector<dcmplx> a54;
        double c2{1.0/2};
};


}