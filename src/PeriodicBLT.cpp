
#include <cmath>
#include "spida/PeriodicBLT.h"
#include "spida/grid/uniformT.h"
#include "spida/transform/fftBLT.h"
#include <pwutils/pwexcept.h>

namespace spida{

  PeriodicBLT::PeriodicBLT(const UniformGridT& grid) :
      m_gr(new UniformGridT(grid)),
      m_tr(new FFTBLT(grid)),
      m_vs(grid.getNst())
  {
  }

  PeriodicBLT::~PeriodicBLT()
  {
      delete m_tr;
      delete m_gr;
  }


  void PeriodicBLT::dT(const std::vector<double>& in,std::vector<double>& out,int n) 
  {
    if(n < 0){
        throw pw::Exception("PeriodicBLT::DT","DT(in,out,n) must have n >= 0 ");
    }
    m_tr->T_To_ST(in,m_vs);
    const std::vector<double>& st = m_gr->getST();
    for(auto i = 0; i < st.size(); i++)
        m_vs[i] = std::pow(ii*st[i],n)*m_vs[i];
    m_tr->ST_To_T(m_vs,out);
  }

}

