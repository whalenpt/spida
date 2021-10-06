
#include <cmath>
#include <memory>
#include "spida/SpidaRCVT.h"
#include "spida/grid/uniformCVT.h"
#include "spida/grid/besselR.h"
#include "spida/transform/hankelfftRCVT.h"
#include <pwutils/pwexcept.h>

namespace spida{

  SpidaRCVT::SpidaRCVT(const BesselRootGridR& gridR,\
          const UniformGridCVT& gridT,unsigned threads) :
      m_gridR(std::make_unique<BesselRootGridR>(gridR)),
      m_gridT(std::make_unique<UniformGridCVT>(gridT)),
      m_tr(std::make_unique<HankelFFTRCVT>(gridR,gridT,threads)) { }
  
  SpidaRCVT::~SpidaRCVT() {}

  const BesselRootGridR& SpidaRCVT::getGridR() const { return *m_gridR; }
  const UniformGridCVT& SpidaRCVT::getGridT() const { return *m_gridT; }
  const HankelFFTRCVT& SpidaRCVT::getTransformRT() const { return *m_tr; }
  const std::vector<double>& SpidaRCVT::getR() const {return m_gridR->getR();}
  const std::vector<double>& SpidaRCVT::getSR() const {return m_gridR->getSR();}
  const std::vector<double>& SpidaRCVT::getT() const {return m_gridT->getT();}
  const std::vector<double>& SpidaRCVT::getST() const {return m_gridT->getST();}
  void SpidaRCVT::RT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->RT_To_SRST(in,out);} 
  void SpidaRCVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->SRST_To_RT(in,out);} 
  unsigned SpidaRCVT::spectralSize() const {return m_gridR->getNr()*m_gridT->getNt();}
  unsigned SpidaRCVT::physicalSize() const {return m_gridR->getNsr()*m_gridT->getNst();}

  void SpidaRCVT::mirrorR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const
  {
      for(auto i = 0; i < m_gridT->getNt(); i++){
          m_gridR->mirrorGrid(in.data()+i*m_gridR->getNr(),out.data()+2*i*m_gridR->getNr());
      }
  }

//  void SpidaRCVT::dT(const std::vector<double>& in,std::vector<double>& out,int n) 
//  {
//    if(n < 0){
//        throw pw::Exception("SpidaRCVT::DT","DT(in,out,n) must have n >= 0 ");
//    }
//    m_tr->T_To_ST(in,m_vs);
//    const std::vector<double>& st = m_gr->getST();
//    for(auto i = 0; i < st.size(); i++)
//        m_vs[i] = std::pow(ii*st[i],n)*m_vs[i];
//    m_tr->ST_To_T(m_vs,out);
//  }

}





