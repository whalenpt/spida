
#include "spida/constants.h"
#include "spida/spidaFFTWX.h"
#include "spida/grid/gridX.h"
#include "spida/grid/chebX.h"
#include "spida/transform/transformX.h"
#include "spida/transform/chebFFTWX.h"
#include <vector>
#include <memory>
#include <iostream>

namespace spida{

  ChebFFTWX::ChebFFTWX(const ChebExtremaGridX& grid) : SpidaX(grid),
      m_gr(new ChebExtremaGridX(grid)),
      m_tr(new ChebTransformFFTWX(grid)),
      m_sp(grid.getNx()), 
      m_dsp(grid.getNx())
  {
  }

  ChebFFTWX::~ChebFFTWX()
  {
      delete m_tr;
      delete m_gr;
  }

  const GridX& ChebFFTWX::getGridX() {
      return *m_gr;
  }
  const TransformX& ChebFFTWX::getTransformX() {
      return *m_tr;
  }

  void ChebFFTWX::X_To_SX(const std::vector<double>& in,std::vector<double>& out) {
      m_tr->X_To_SX(in,out);
  } 
  void ChebFFTWX::SX_To_X(const std::vector<double>& in,std::vector<double>& out) {
      m_tr->SX_To_X(in,out);
  } 

  void ChebFFTWX::dX(const std::vector<double>& in,std::vector<double>& out,int n) 
  {
      X_To_SX(in,m_sp);
      dSX(m_sp,n);
      SX_To_X(m_sp,out);
  }

  void ChebFFTWX::dX(const std::vector<double>& in,std::vector<double>& out) 
  {
      ChebFFTWX::dX(in,out,1);
  }

  void ChebFFTWX::dSX(std::vector<double>& a,int n)
  {
      if(n < 0){ 
          throw std::invalid_argument("ChebFFTWX::dSX(std::vector<double>& a,int n) \
                  n must be a non-negative integer");
      } else if(n > 8){
          throw std::invalid_argument("ChebFFTWX::dSX(std::vector<double>& a,int n) \
                  n must be less than or equal to 8");
      }
      for(int i = 0; i < n; i++)
          ChebFFTWX::dSX(a);
  }

  void ChebFFTWX::dSX(std::vector<double>& a)
  {
      int N = m_gr->getNsx();
      // Scaling factor for grid size (wider grid -> narrower spectrum and vice versa)
      double sf = 2.0/m_gr->getL();
      m_dsp[N-1] = 0.0;
      m_dsp[N-2] =  2.0*(N-1)*sf*a[N-1]; 
      for(int k = N-2; k > 1; k--)
          m_dsp[k-1] = m_dsp[k+1] + 2.0*sf*k*a[k];
      m_dsp[0] = 0.5*m_dsp[2] + sf*a[1];
      std::copy(std::begin(m_dsp),std::end(m_dsp),std::begin(a));
  }


  void ChebFFTWX::dSX(const std::vector<double>& a,std::vector<double>& b,int n) 
  {
      std::copy(std::cbegin(a),std::cend(a),std::begin(b));
      dSX(b,n);
  }

  void ChebFFTWX::dSX(const std::vector<double>& a,std::vector<double>& b) 
  {
      std::copy(std::cbegin(a),std::cend(a),std::begin(b));
      dSX(b);
  }

}




