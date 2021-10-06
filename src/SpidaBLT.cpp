
#include <cmath>
#include <memory>
#include "spida/SpidaBLT.h"
#include "spida/grid/uniformT.h"
#include "spida/transform/fftBLT.h"
#include <pwutils/pwexcept.h>

namespace spida{

  SpidaBLT::SpidaBLT(const UniformGridRVT& grid) :
      m_gr(std::make_unique<UniformGridRVT>(grid)),
      m_tr(std::make_unique<FFTBLT>(grid)),
      m_vs(grid.getNst()) {}

  SpidaBLT::~SpidaBLT() {}
  const UniformGridRVT& SpidaBLT::getGridT() const { return *m_gr; }
  const FFTBLT& SpidaBLT::getTransformT() const { return *m_tr; }
  const std::vector<double>& SpidaBLT::getT() const {return m_gr->getT();}
  const std::vector<double>& SpidaBLT::getST() const {return m_gr->getST();}
  void SpidaBLT::T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out) {m_tr->T_To_ST(in,out);} 
  void SpidaBLT::ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out) {m_tr->ST_To_T(in,out);} 
  void SpidaBLT::CVT_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->CVT_To_ST(in,out);} 
  void SpidaBLT::ST_To_CVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->ST_To_CVT(in,out);} 

  void SpidaBLT::dT(const std::vector<double>& in,std::vector<double>& out,int n) 
  {
    if(n < 0){
        throw pw::Exception("SpidaBLT::DT","DT(in,out,n) must have n >= 0 ");
    }
    m_tr->T_To_ST(in,m_vs);
    const std::vector<double>& st = m_gr->getST();
    for(auto i = 0; i < st.size(); i++)
        m_vs[i] = std::pow(ii*st[i],n)*m_vs[i];
    m_tr->ST_To_T(m_vs,out);
  }

}

