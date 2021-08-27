

#include "spida/solver/solverETDCS.h"
#include "spida/model/model.h"
#include "spida/propagator/propagator.h"
#include <cmath>

namespace spida{

SolverCS_ETD::SolverCS_ETD(ModelCV& cmodel,PropagatorCV& cprs)
 :  SolverCV_CS(cmodel,cprs)
{
    m_mode_cutoff = 0.01;
    m_contour_radi = 1.0;
    m_contourM = 32;
}

ETD4::ETD4(ModelCV& model,PropagatorCV& prs) :
  SolverCS_ETD(model,prs),sz(SolverCV::size()),
    L(model.linOp()), EL(sz), EL2(sz), N1(sz), N2(sz),
    N3(sz), N4(sz), tempK(sz), a21(sz), a31(sz), a32(sz), a41(sz), a43(sz),
    a51(sz), a52(sz), a54(sz)
{
    c1 = 0.0; c2 = 1.0/2; c3 = 1.0/2; c4 = 1.0;
}

void ETD4::updateCoefficients(double dt)
{
    std::vector<dcmplx> r(contourPoints(),0.0);
    for(int j = 0; j < contourPoints(); j++){
        double expv = 2.0*PI*(j+0.5)/contourPoints();
        r[j] = contourRadius()*exp(ii*expv);
    }

    for(int i = 0; i < sz; i++) {
        dcmplx Lval(dt*L[i]);
        EL[i] = exp(Lval);
        EL2[i] = exp(c2*Lval);
        dcmplx a21val(0.0,0.0); dcmplx a31val(0.0,0.0); dcmplx a32val(0.0,0.0);
        dcmplx a41val(0.0,0.0); dcmplx a43val(0.0,0.0); 
        dcmplx a51val(0.0,0.0); dcmplx a52val(0.0,0.0); dcmplx a54val(0.0,0.0);
        if(abs(Lval) < modeCutoff()){
            for(int j = 0; j < contourPoints(); j++){
                dcmplx z(Lval + r[j]);
                a21val += (-1.0 + exp(z/2.0))/z;
                a31val += (4.0 + exp(z/2.0)*(-4.0 + z) + z)/pow(z,2);
                a32val += (-2.0*(2.0 - 2.0*exp(z/2.0) + z))/pow(z,2);
                a41val += (2.0 + exp(z)*(-2.0 + z) + z)/pow(z,2);
                a43val += (-2.0*(1.0 - exp(z) + z))/pow(z,2);
                a51val += (-4.0 - z + exp(z)*(4.0 - 3.0*z + pow(z,2)))/pow(z,3);
                a52val += (2.0*(2.0 + exp(z)*(-2.0 + z) + z))/pow(z,3);
                a54val += -((4.0 + exp(z)*(-4.0 + z) + 3.0*z + pow(z,2))/pow(z,3));
            }
            double mult = dt/contourPoints();
            a21[i] = mult*a21val; a31[i] = mult*a31val; a32[i] = mult*a32val;
            a41[i] = mult*a41val; a43[i] = mult*a43val;
            a51[i] = mult*a51val; a52[i] = mult*a52val; a54[i] = mult*a54val;
        }
        else{
            a21[i] = dt*(-1.0 + exp(Lval/2.0))/Lval;
            a31[i] = dt*(4.0 + exp(Lval/2.0)*(-4.0 + Lval) + Lval)/pow(Lval,2);
            a32[i] = dt*(-2.0*(2.0 - 2.0*exp(Lval/2.0) + Lval))/pow(Lval,2);
            a41[i] = dt*(2.0 + exp(Lval)*(-2.0 + Lval) + Lval)/pow(Lval,2);
            a43[i] = dt*(-2.0*(1.0 - exp(Lval) + Lval))/pow(Lval,2);
            a51[i] = dt*(-4.0 - Lval + exp(Lval)*(4.0 - 3.0*Lval + pow(Lval,2)))/pow(Lval,3);
            a52[i] = dt*(2.0*(2.0 + exp(Lval)*(-2.0 + Lval) + Lval))/pow(Lval,3);
            a54[i] = dt*-((4.0 + exp(Lval)*(-4.0 + Lval) + 3.0*Lval + pow(Lval,2))/pow(Lval,3));
        }
    }
}

void ETD4::updateStages(std::vector<dcmplx>& in)
{
    SolverCV::model().nonLinResponse(in,N1);
    for(int i = 0; i < sz; i++)
        tempK[i] = EL2[i]*in[i] + a21[i]*N1[i];
  
    SolverCV::model().nonLinResponse(tempK,N2);
    for(int i = 0; i < sz; i++)
        tempK[i] = EL2[i]*in[i] + a31[i]*N1[i] + a32[i]*N2[i];

    SolverCV::model().nonLinResponse(tempK,N3);
    for(int i = 0; i < sz; i++)
        tempK[i] = EL[i]*in[i] + a41[i]*N1[i] + a43[i]*N3[i];

    SolverCV::model().nonLinResponse(tempK,N4);
    for(int i = 0; i < sz; i++) 
        in[i] = EL[i]*in[i] + a51[i]*N1[i] + a52[i]*(N2[i]+N3[i]) + a54[i]*N4[i];
}

}


