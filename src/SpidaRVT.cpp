
#include <cmath>
#include <memory>
#include "spida/RVT.h"
#include "spida/grid/uniformRVT.h"
#include "spida/transform/fftRVT.h"
#include <pwutils/pwexcept.h>

namespace spida{

  SpidaRVT::SpidaRVT(const UniformGridRVT& grid) :
      m_gr(std::make_unique<UniformGridRVT>(grid)),
      m_tr(std::make_unique<FFTRVT>(grid)),
      m_vs(grid.getNst()) {}

  const UniformGridRVT& SpidaRVT::getGridT() const { return *m_gr; }
  const FFTRVT& SpidaRVT::getTransformT() const { return *m_tr; }
  const std::vector<double>& SpidaRVT::getT() const {return m_gr->getT();}
  const std::vector<double>& SpidaRVT::getST() const {return m_gr->getST();}
  void SpidaRVT::T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out) {m_tr->T_To_ST(in,out);} 
  void SpidaRVT::ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out) {m_tr->ST_To_T(in,out);} 
  void SpidaRVT::CVT_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->CVT_To_ST(in,out);} 
  void SpidaRVT::ST_To_CVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->ST_To_CVT(in,out);} 

  void SpidaRVT::dT(const std::vector<double>& in,std::vector<double>& out,int n) 
  {
    if(n < 0){
        throw pw::Exception("SpidaRVT::DT","DT(in,out,n) must have n >= 0 ");
    }
    m_tr->T_To_ST(in,m_vs);
    const std::vector<double>& st = m_gr->getST();
    for(auto i = 0; i < st.size(); i++)
        m_vs[i] = std::pow(-ii*st[i],n)*m_vs[i];
    m_tr->ST_To_T(m_vs,out);
  }

}

