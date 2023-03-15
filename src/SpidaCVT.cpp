#include <cmath>
#include <memory>
#include "spida/CVT.h"
#include "spida/grid/uniformCVT.h"
#include "spida/transform/fftCVT.h"

namespace spida{

  SpidaCVT::SpidaCVT(const UniformGridCVT& grid) :
      m_gr(std::make_unique<UniformGridCVT>(grid)),
      m_tr(std::make_unique<FFTCVT>(grid)),
      m_vs(grid.getNst()) {}

  const UniformGridCVT& SpidaCVT::getGridT() const { return *m_gr; }
  const FFTCVT& SpidaCVT::getTransformT() const { return *m_tr; }
  const std::vector<double>& SpidaCVT::getT() const {return m_gr->getT();}
  const std::vector<double>& SpidaCVT::getST() const {return m_gr->getST();}
  void SpidaCVT::T_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {m_tr->T_To_ST(in,out);} 
  void SpidaCVT::ST_To_T(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {m_tr->ST_To_T(in,out);} 

  void SpidaCVT::dT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n) noexcept
  {
    if(n == 0)
        return;
    m_tr->T_To_ST(in,out);
    dST(out,m_vs,n);
    m_tr->ST_To_T(m_vs,out);
  }

  void SpidaCVT::dST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n) noexcept
  {
    if(n == 0)
        return;
    const std::vector<double>& omega = m_gr->getST();
    for(size_t i = 0; i < omega.size(); i++)
        out[i] = std::pow(-ii*omega[i],n)*in[i];
  }

}






