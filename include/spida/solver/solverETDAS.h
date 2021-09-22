#ifndef SOLVERETDAS_H_
#define SOLVERETDAS_H_

#include "spida/model/model.h"
#include "spida/solver/solver.h"
#include "spida/propagator/propagator.h"
#include "spida/helper/constants.h"
#include <iostream>

//class Model;
namespace spida{

class SolverAS_ETD : public SolverCV_AS
{
    public:
        SolverAS_ETD(ModelCV* cmodel,double sf,double qv);
        virtual ~SolverAS_ETD() {}; 
        void setModeCutoff(double val) {m_mode_cutoff = val;}
        void setContourRadius(double val) {m_contour_radi = val;}
        void setContourPoints(int val) {m_contourM = val;}
        double modeCutoff() {return m_mode_cutoff;}
        double contourRadius() {return m_contour_radi;}
        int contourPoints() {return m_contourM;}
    private:
        virtual void updateCoefficients(double dt) noexcept = 0;
        virtual void updateStages(const std::vector<dcmplx>& in,\
                std::vector<dcmplx>& y,std::vector<dcmplx>& err) noexcept = 0;
        double m_mode_cutoff;
        int m_contourM;
        double m_contour_radi;
};

class ETD34 : public SolverAS_ETD
{
    public:
        ETD34(ModelCV* cmodel);
        ~ETD34() {};
    private:
        void updateCoefficients(double dt) noexcept;
        void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,\
                 std::vector<dcmplx>& err) noexcept;
        int m_sz;
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
        double c1,c2,c3,c4,c5;
        bool N1_init;
        pw::ThreadManager& thmgt;
        void worker_coeff(double ds,int tid);
};


class ETD35: public SolverAS_ETD
{
    public:
        ETD35(ModelCV* cmodel);
        ~ETD35() {};
    private:
        void updateCoefficients(double dt) noexcept;
        void updateStages(const std::vector<dcmplx>& in,\
                std::vector<dcmplx>& y,std::vector<dcmplx>& err) noexcept;

        int m_sz;
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
        double c1,c2,c3,c4,c5,c6,c7;
        bool N1_init;
        pw::ThreadManager& thmgt;
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

#endif




