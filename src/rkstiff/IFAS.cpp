
#include "spida/rkstiff/IFAS.h"
#include <cmath>

namespace spida{

IF34::IF34(const LinOp& Lop,const NLfunc& NL,double sf,double qv,bool use_refs)
 :  SolverCV_AS(Lop,NL,sf,qv,use_refs), m_sz(SolverCV::size()),
    L(SolverCV::L()),EL(m_sz), EL2(m_sz), 
    N1(m_sz), N2(m_sz), N3(m_sz), N4(m_sz), N5(m_sz), 
    tempK(m_sz), a21(m_sz), a43(m_sz), 
    a51(m_sz), a52(m_sz)
{
  statCenter().setHeader("IF34 STATS");
  statCenter().addCounter("Nonlinear Function Evaluations");
}

void IF34::updateCoefficients(double dt) noexcept
{
    for(unsigned i = 0; i < m_sz; i++) {
        EL[i] = exp(dt*L[i]);
        EL2[i] = exp(c2*dt*L[i]);
        a21[i] = (1.0/2.0)*dt*exp(dt*L[i]/2.0);
        a43[i] = dt*exp(dt*L[i]/2.0);
        a51[i] = (1.0/6.0)*dt*exp(dt*L[i]);
        a52[i] = (1.0/3.0)*dt*exp(dt*L[i]/2.0);
    }
    a32 = (1.0/2.0)*dt;
    a54 = (1.0/6.0)*dt;
}

void IF34::updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& ynew,\
        std::vector<dcmplx>& errVec) noexcept
{
    if(!N1_init){
        SolverCV::NL()(in,N1);
        N1_init = true;
        statCenter().incrementCounter("Nonlinear Function Evaluations",1);
    }
    else if(SolverCV_AS::accept()){
        for(unsigned i = 0; i < m_sz; i++) 
            N1[i] = N5[i]; 
    }
    for(unsigned i = 0; i < m_sz; i++)
        tempK[i] =  EL2[i]*in[i] + a21[i]*N1[i];
    SolverCV::NL()(tempK,N2);

    for(unsigned i = 0; i < m_sz; i++) 
        tempK[i] = EL2[i]*in[i] + a32*N2[i];
    SolverCV::NL()(tempK,N3);

    for(unsigned i = 0; i < m_sz; i++) 
        tempK[i] = EL[i]*in[i] + a43[i]*N3[i];
    SolverCV::NL()(tempK,N4);

    for(unsigned i = 0; i < m_sz; i++) 
        ynew[i] = EL[i]*in[i] + a51[i]*N1[i] + a52[i]*(N2[i]+N3[i]) + a54*N4[i];
    SolverCV::NL()(ynew,N5);

    for (unsigned i = 0; i < m_sz; i++)
        errVec[i] = a54*(N4[i] - N5[i]);

    statCenter().incrementCounter("Nonlinear Function Evaluations",4);
}

IF45DP::IF45DP(const LinOp& Lop,const NLfunc& NL,double sf,double qv,bool use_refs)
 :  SolverCV_AS(Lop,NL,sf,qv,use_refs), m_sz(SolverCV::size()),
    L(SolverCV::L()),EL(m_sz), EL2(m_sz), EL3(m_sz), EL4(m_sz), EL5(m_sz),
    N1(m_sz), N2(m_sz), N3(m_sz), N4(m_sz), N5(m_sz), N6(m_sz), N7(m_sz),
    tempK(m_sz), a21(m_sz), a31(m_sz), a32(m_sz), a41(m_sz), a42(m_sz), a43(m_sz), 
    a51(m_sz), a52(m_sz), a53(m_sz), a54(m_sz),
    a61(m_sz), a62(m_sz), a63(m_sz), a64(m_sz), a65(m_sz),
    a71(m_sz), a73(m_sz), a74(m_sz), a75(m_sz),
    r1(m_sz), r3(m_sz), r4(m_sz), r5(m_sz)
{
    SolverCV::statCenter().setHeader("IF45DP STATS");
    SolverCV::statCenter().addCounter("Nonlinear Function Evaluations");
}

void IF45DP::updateCoefficients(double dt) noexcept
{
    for(unsigned i = 0; i < m_sz; i++) {
        EL[i] = exp(dt*L[i]);
        EL2[i] = exp(c2*dt*L[i]);
        EL3[i] = exp(c3*dt*L[i]);
        EL4[i] = exp(c4*dt*L[i]);
        EL5[i] = exp(c5*dt*L[i]);
        a21[i] = (1.0/5.0)*dt*exp(c2*dt*L[i]);
        a31[i] = (3.0/40.0)*dt*exp(c3*dt*L[i]);
        a32[i] = (9.0/40.0)*dt*exp(dt*L[i]/10.0);
        a41[i] = (44.0/45.0)*dt*exp(c4*dt*L[i]);
        a42[i] = (-56.0/15.0)*dt*exp(3.0*dt*L[i]/5.0);
        a43[i] = (32.0/9.0)*dt*exp(dt*L[i]/2.0);
        a51[i] = (19372.0/6561.0)*dt*exp(c5*dt*L[i]);
        a52[i] = (-25360.0/2187.0)*dt*exp(31.0*dt*L[i]/45.0);
        a53[i] = (64448.0/6561.0)*dt*exp(53.0*dt*L[i]/90.0);
        a54[i] = (-212.0/729.0)*dt*exp(4.0*dt*L[i]/45.0);
        a61[i] = (9017.0/3168.0)*dt*exp(dt*L[i]);
        a62[i] = (-355.0/33.0)*dt*exp(4.0*dt*L[i]/5.0);
        a63[i] = (46732.0/5247.0)*dt*exp(7.0*dt*L[i]/10.0);
        a64[i] = (49.0/176.0)*dt*exp(dt*L[i]/5.0);
        a65[i] = (-5103.0/18656.0)*dt*exp(dt*L[i]/9.0);
        a71[i] = (35.0/384.0)*dt*exp(dt*L[i]);
        a73[i] = (500.0/1113.0)*dt*exp(7.0*dt*L[i]/10.0);
        a74[i] = (125.0/192.0)*dt*exp(dt*L[i]/5.0);
        a75[i] = (-2187.0/6784.0)*dt*exp(dt*L[i]/9.0);
        r1[i] = dt*71.0*exp(dt*L[i])/57600.0;
        r3[i] = -dt*71.0*exp(7.0*dt*L[i]/10.0)/16695.0;
        r4[i] = dt*17.0*exp(dt*L[i]/5.0)/1920.0;
        r5[i] = -dt*17253.0*exp(dt*L[i]/9.0)/339200.0;
    }
    a76 = (11.0/84.0)*dt;
    r6 = (22.0/525.0)*dt;
    r7 = (-1.0/40.0)*dt;
}

void IF45DP::updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& ynew,\
        std::vector<dcmplx>& errVec) noexcept
{
    if(!N1_init){
        SolverCV::NL()(in,N1);
        N1_init = true;
        statCenter().incrementCounter("Nonlinear Function Evaluations",1);
    }
    else if(SolverCV_AS::accept()){
        for(unsigned i = 0; i < m_sz; i++) 
            N1[i] = N7[i]; 
    }
    for(unsigned i = 0; i < m_sz; i++)
        tempK[i] =  EL2[i]*in[i] + a21[i]*N1[i];
    SolverCV::NL()(tempK,N2);

    for(unsigned i = 0; i < m_sz; i++) 
        tempK[i] = EL3[i]*in[i] + a31[i]*N1[i] + a32[i]*N2[i];
    SolverCV::NL()(tempK,N3);

    for(unsigned i = 0; i < m_sz; i++) 
        tempK[i] = EL4[i]*in[i] + a41[i]*N1[i] + a42[i]*N2[i] + a43[i]*N3[i];
    SolverCV::NL()(tempK,N4);

    for(unsigned i = 0; i < m_sz; i++) 
        tempK[i] = EL5[i]*in[i] + a51[i]*N1[i] + a52[i]*N2[i] + a53[i]*N3[i]+ a54[i]*N4[i];
    SolverCV::NL()(tempK,N5);

    for(unsigned i = 0; i < m_sz; i++) 
        tempK[i] = EL[i]*in[i] + a61[i]*N1[i] + a62[i]*N2[i] + a63[i]*N3[i] \
                   + a64[i]*N4[i] + a65[i]*N5[i];
    SolverCV::NL()(tempK,N6);

    for(unsigned i = 0; i < m_sz; i++) 
        ynew[i] = EL[i]*in[i] + a71[i]*N1[i] + a73[i]*N3[i] \
                   + a74[i]*N4[i] + a75[i]*N5[i] + a76*N6[i];
    SolverCV::NL()(tempK,N7);

    for (unsigned i = 0; i < m_sz; i++)
        errVec[i] = r1[i]*N1[i]+r3[i]*N3[i] + r4[i]*N4[i] + r5[i]*N5[i] + r6*N6[i] + r7*N7[i];

    statCenter().incrementCounter("Nonlinear Function Evaluations",6);
}

}