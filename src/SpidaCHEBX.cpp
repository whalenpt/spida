

#include "spida/chebX.h"
#include "spida/grid/chebX.h"
#include "spida/transform/chebX.h"
#include "spida/helper/constants.h"
#include <vector>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <cassert>


namespace spida{

  SpidaChebX::SpidaChebX(const ChebRootGridX& grid) :
    m_gr(std::make_unique<ChebRootGridX>(grid)),
    m_tr(std::make_unique<ChebTransformX>(grid)),
    m_sp(grid.getNx()), 
    m_dsp(grid.getNsx()) { }

  const GridX& SpidaChebX::getGridX() {
      return *m_gr;
  }
  const ChebTransformX& SpidaChebX::getTransformX() {
      return *m_tr;
  }

  void SpidaChebX::X_To_SX(const std::vector<double>& in,std::vector<double>& out) {
      m_tr->X_To_SX(in,out);
  } 
  void SpidaChebX::SX_To_X(const std::vector<double>& in,std::vector<double>& out) {
      m_tr->SX_To_X(in,out);
  } 

  void SpidaChebX::dX(const std::vector<double>& in,std::vector<double>& out,int n) 
  {
      X_To_SX(in,m_sp);
      dSX(m_sp,n);
      SX_To_X(m_sp,out);
  }

  void SpidaChebX::dX(const std::vector<double>& in,std::vector<double>& out) 
  {
      SpidaChebX::dX(in,out,1);
  }

  void SpidaChebX::dSX(std::vector<double>& a,int n)
  {
      if(n < 0){ 
          throw std::invalid_argument("SpidaChebX::dSX(std::vector<double>& a,int n) "
                  "n must be a non-negative integer");
      } else if(n > 8){
          throw std::invalid_argument("SpidaChebX::dSX(std::vector<double>& a,int n) "
                  "n must be less than or equal to 8");
      }
      for(int i = 0; i < n; i++)
          SpidaChebX::dSX(a);
  }

  void SpidaChebX::dSX(std::vector<double>& a)
  {
      int N = m_gr->getNsx();
      double sf = 2.0/m_gr->getL();
      m_dsp[N-1] = 0.0;
      m_dsp[N-2] =  2.0*(N-1)*sf*a[N-1]; 
      for(int k = N-2; k > 0; k--)
          m_dsp[k-1] = m_dsp[k+1] + 2.0*sf*k*a[k];
      std::copy(std::begin(m_dsp),std::end(m_dsp),std::begin(a));
  }

  void SpidaChebX::dSX(const std::vector<double>& a,std::vector<double>& b,int n) 
  {
      std::copy(std::cbegin(a),std::cend(a),std::begin(b));
      dSX(b,n);
  }

  void SpidaChebX::dSX(const std::vector<double>& a,std::vector<double>& b) 
  {
      std::copy(std::cbegin(a),std::cend(a),std::begin(b));
      dSX(b);
  }

}