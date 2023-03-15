// ETDAS.h
#pragma once

#include <iostream>
#include "spida/rkstiff/solver.h"
#include "spida/helper/constants.h"

namespace spida{

class SolverAS_ETD : public SolverCV_AS
{
    public:
        using SolverCV_AS::SolverCV_AS;
        ~SolverAS_ETD() override = default;
        void setModeCutoff(double val) {m_mode_cutoff = val;}
        void setContourRadius(double val) {m_contour_radi = val;}
        void setContourPoints(unsigned val) {m_contourM = val;}
        double modeCutoff() const {return m_mode_cutoff;}
        double contourRadius() const {return m_contour_radi;}
        unsigned contourPoints() const {return m_contourM;}
    private:
        double m_mode_cutoff{0.01};
        unsigned m_contourM{32};
        double m_contour_radi{1.0};
};

class ETD34 : public SolverAS_ETD
{
    public:
        ETD34(const LinOp& Lop,const NLfunc& NL,bool use_refs=false);
        ~ETD34() override = default;
    private:
        void updateCoefficients(double dt) noexcept override;
        void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,\
                 std::vector<dcmplx>& err) noexcept override;
        unsigned m_sz;
        const std::vector<dcmplx>& L;
        std::vector<dcmplx> EL; std::vector<dcmplx> EL2; 
        std::vector<dcmplx> N1; std::vector<dcmplx> N2; std::vector<dcmplx> N3; 
        std::vector<dcmplx> N4; std::vector<dcmplx> N5;
        std::vector<dcmplx> tempK;
        std::vector<dcmplx> a21;
        std::vector<dcmplx> a31; std::vector<dcmplx> a32;
        std::vector<dcmplx> a41; std::vector<dcmplx> a43;
        std::vector<dcmplx> a51; std::vector<dcmplx> a52; std::vector<dcmplx> a54;
        std::vector<dcmplx> r1; 
        double c2{1.0/2};
        bool N1_init{false};
        void worker_coeff(double ds,int tid);
};


class ETD35: public SolverAS_ETD
{
    public:
        ETD35(const LinOp& Lop,const NLfunc& NL,bool use_refs = false);
        ~ETD35() override = default;
    private:
        void updateCoefficients(double dt) noexcept override;
        void updateStages(const std::vector<dcmplx>& in,\
                std::vector<dcmplx>& y,std::vector<dcmplx>& err) noexcept override;

        unsigned m_sz;
        const std::vector<dcmplx>& L;
        std::vector<dcmplx> EL; std::vector<dcmplx> EL2; std::vector<dcmplx> EL4;
        std::vector<dcmplx> EL5; 
        std::vector<dcmplx> N1; std::vector<dcmplx> N2; std::vector<dcmplx> N3; std::vector<dcmplx> N4;
        std::vector<dcmplx> N5; std::vector<dcmplx> N6; std::vector<dcmplx> N7;
        std::vector<dcmplx> tempK;
        std::vector<dcmplx> a21;
        std::vector<dcmplx> a31; std::vector<dcmplx> a32;
        std::vector<dcmplx> a41; std::vector<dcmplx> a43;
        std::vector<dcmplx> a51; std::vector<dcmplx> a52; std::vector<dcmplx> a54;
        std::vector<dcmplx> a61; std::vector<dcmplx> a62; std::vector<dcmplx> a63; std::vector<dcmplx> a65;
        std::vector<dcmplx> a71; std::vector<dcmplx> a73; std::vector<dcmplx> a74; std::vector<dcmplx> a75; std::vector<dcmplx> a76;

        double c2{1.0/4}; double c4{1.0/2}; double c5{3.0/4};
        bool N1_init{false};
        void worker_coeff(double ds,int tid);
        void worker_stage2(const std::vector<dcmplx>& in,int sti,int endi);
        void worker_stage3(const std::vector<dcmplx>& in,int sti,int endi);
        void worker_stage4(const std::vector<dcmplx>& in,int sti,int endi);
        void worker_stage5(const std::vector<dcmplx>& in,int sti,int endi);
        void worker_stage6(const std::vector<dcmplx>& in,int sti,int endi);
        void worker_stage7(const std::vector<dcmplx>& in,std::vector<dcmplx>& ynew,\
        std::vector<dcmplx>& errVec,int sti,int endi);

};

}