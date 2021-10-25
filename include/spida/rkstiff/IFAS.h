
#pragma once

#include "spida/rkstiff/solver.h"
#include "spida/helper/constants.h"

namespace spida{

class IF34 : public SolverCV_AS
{
  public:
    IF34(const LinOp& Lop,const NLfunc& NL,double sf = 0.84,double qv = 4.0);
    ~IF34() {};
  private:
    void updateCoefficients(double dt) noexcept;
    void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,\
             std::vector<dcmplx>& err) noexcept;

    unsigned m_sz;
    const std::vector<dcmplx>& L;
    std::vector<dcmplx> EL; std::vector<dcmplx> EL2; 
    std::vector<dcmplx> N1; std::vector<dcmplx> N2; std::vector<dcmplx> N3; 
    std::vector<dcmplx> N4; std::vector<dcmplx> N5;
    std::vector<dcmplx> tempK;

    double a32; double a54;
    std::vector<dcmplx> a21;
    std::vector<dcmplx> a43;
    std::vector<dcmplx> a51; std::vector<dcmplx> a52; 

    double c1,c2,c3,c4,c5; //c6 = 1.0;
    bool N1_init;
};

class IF45DP : public SolverCV_AS
{
  public:
    IF45DP(const LinOp& Lop,const NLfunc& NL,double sf = 0.9,double qv = 5.0);
    ~IF45DP() {};
  private:
    void updateCoefficients(double dt) noexcept;
    void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,\
             std::vector<dcmplx>& err) noexcept;

    unsigned m_sz;
    const std::vector<dcmplx>& L;
    std::vector<dcmplx> EL; std::vector<dcmplx> EL2; 
    std::vector<dcmplx> EL3; std::vector<dcmplx> EL4; std::vector<dcmplx> EL5; 
    std::vector<dcmplx> N1; std::vector<dcmplx> N2; std::vector<dcmplx> N3; 
    std::vector<dcmplx> N4; std::vector<dcmplx> N5;
    std::vector<dcmplx> N6; std::vector<dcmplx> N7;
    std::vector<dcmplx> tempK;

    double a76; double r6; double r7;
    std::vector<dcmplx> a21;
    std::vector<dcmplx> a31; std::vector<dcmplx> a32;
    std::vector<dcmplx> a41; std::vector<dcmplx> a42; std::vector<dcmplx> a43;
    std::vector<dcmplx> a51; std::vector<dcmplx> a52; std::vector<dcmplx> a53; std::vector<dcmplx> a54;
    std::vector<dcmplx> a61; std::vector<dcmplx> a62; std::vector<dcmplx> a63; std::vector<dcmplx> a64; std::vector<dcmplx> a65;
    std::vector<dcmplx> a71; std::vector<dcmplx> a73; std::vector<dcmplx> a74; std::vector<dcmplx> a75; 
    std::vector<dcmplx> r1; std::vector<dcmplx> r3; std::vector<dcmplx> r4; std::vector<dcmplx> r5;
    double c1,c2,c3,c4,c5,c6,c7;
    bool N1_init;
};




}



