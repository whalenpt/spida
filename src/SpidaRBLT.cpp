
#include <cmath>
#include <memory>
#include "spida/SpidaRBLT.h"
#include "spida/grid/uniformT.h"
#include "spida/grid/besselR.h"
#include "spida/transform/hankelfftRBLT.h"
#include <pwutils/pwexcept.h>

namespace spida{

  SpidaRBLT::SpidaRBLT(const BesselRootGridR& gridR,const UniformGridT& gridT,unsigned int threads) :
      m_gridR(std::make_unique<BesselRootGridR>(gridR)),
      m_gridT(std::make_unique<UniformGridT>(gridT)),
      m_tr(std::make_unique<HankelFFTRBLT>(gridR,gridT,threads)) { }
  
  SpidaRBLT::~SpidaRBLT() {}

  const BesselRootGridR& SpidaRBLT::getGridR() const { return *m_gridR; }
  const UniformGridT& SpidaRBLT::getGridT() const { return *m_gridT; }
  const HankelFFTRBLT& SpidaRBLT::getTransformRT() const { return *m_tr; }
  const std::vector<double>& SpidaRBLT::getR() const {return m_gridR->getR();}
  const std::vector<double>& SpidaRBLT::getSR() const {return m_gridR->getSR();}
  const std::vector<double>& SpidaRBLT::getT() const {return m_gridT->getT();}
  const std::vector<double>& SpidaRBLT::getST() const {return m_gridT->getST();}
  void SpidaRBLT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) {m_tr->RT_To_SRST(in,out);} 
  void SpidaRBLT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) {m_tr->SRST_To_RT(in,out);} 


//  void SpidaRBLT::dT(const std::vector<double>& in,std::vector<double>& out,int n) 
//  {
//    if(n < 0){
//        throw pw::Exception("SpidaRBLT::DT","DT(in,out,n) must have n >= 0 ");
//    }
//    m_tr->T_To_ST(in,m_vs);
//    const std::vector<double>& st = m_gr->getST();
//    for(auto i = 0; i < st.size(); i++)
//        m_vs[i] = std::pow(ii*st[i],n)*m_vs[i];
//    m_tr->ST_To_T(m_vs,out);
//  }

}





