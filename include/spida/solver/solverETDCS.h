
#ifndef SOLVERETDCS_H_
#define SOLVERETDCS_H_

#include "spida/solver/solver.h"
#include "spida/helper/constants.h"

namespace spida{

class ModelCV;
class PropagatorCV;

class SolverCS_ETD : public SolverCV_CS
{
    public:
        SolverCS_ETD(ModelCV& cmodel);
        virtual ~SolverCS_ETD() {}; 
        void setModeCutoff(double val) {m_mode_cutoff = val;}
        void setContourRadius(double val) {m_contour_radi = val;}
        void setContourPoints(int val) {m_contourM = val;}
        double modeCutoff() {return m_mode_cutoff;}
        double contourRadius() {return m_contour_radi;}
        int contourPoints() {return m_contourM;}
    private:
        virtual void updateCoefficients(double dt) = 0;
        virtual void updateStages(std::vector<dcmplx>& in) = 0;
        double m_mode_cutoff;
        int m_contourM;
        double m_contour_radi;
};

class ETD4 : public SolverCS_ETD
{
    public:
        ETD4(ModelCV& cmodel);
        ~ETD4() {};
    private:
        void updateCoefficients(double dt);
        void updateStages(std::vector<dcmplx>& in);
        int sz;
        const std::vector<dcmplx>& L;
        std::vector<dcmplx> EL; std::vector<dcmplx> EL2; 
        std::vector<dcmplx> N1; std::vector<dcmplx> N2; 
        std::vector<dcmplx> N3; std::vector<dcmplx> N4;
        std::vector<dcmplx> tempK;
        std::vector<dcmplx> a21; std::vector<dcmplx> a31; std::vector<dcmplx> a32;
        std::vector<dcmplx> a41; std::vector<dcmplx> a43;
        std::vector<dcmplx> a51; std::vector<dcmplx> a52; std::vector<dcmplx> a54;
        double c1,c2,c3,c4;
};


}

#endif


