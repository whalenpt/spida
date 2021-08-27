

#include "spida/ChebX.h"
#include "spida/grid/chebX.h"
#include "spida/transform/chebX.h"
#include "spida/helper/constants.h"
#include <vector>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <cassert>


namespace spida{

  ChebX::ChebX(const ChebRootGridX& grid) :
    m_gr(new ChebRootGridX(grid)),
    m_tr(new ChebTransformX(grid)),
    m_sp(grid.getNx()), 
    m_dsp(grid.getNsx())
  {
  }

  ChebX::~ChebX()
  {
      delete m_tr;
      delete m_gr;
  }

  const GridX& ChebX::getGridX() {
      return *m_gr;
  }
  const ChebTransformX& ChebX::getTransformX() {
      return *m_tr;
  }

  void ChebX::X_To_SX(const std::vector<double>& in,std::vector<double>& out) {
      m_tr->X_To_SX(in,out);
  } 
  void ChebX::SX_To_X(const std::vector<double>& in,std::vector<double>& out) {
      m_tr->SX_To_X(in,out);
  } 

  void ChebX::dX(const std::vector<double>& in,std::vector<double>& out,int n) 
  {
      X_To_SX(in,m_sp);
      dSX(m_sp,n);
      SX_To_X(m_sp,out);
  }

  void ChebX::dX(const std::vector<double>& in,std::vector<double>& out) 
  {
      ChebX::dX(in,out,1);
  }

  void ChebX::dSX(std::vector<double>& a,int n)
  {
      if(n < 0){ 
          throw std::invalid_argument("ChebX::dSX(std::vector<double>& a,int n) \
                  n must be a non-negative integer");
      } else if(n > 8){
          throw std::invalid_argument("ChebX::dSX(std::vector<double>& a,int n) \
                  n must be less than or equal to 8");
      }
      for(int i = 0; i < n; i++)
          ChebX::dSX(a);
  }

  void ChebX::dSX(std::vector<double>& a)
  {
      int N = m_gr->getNsx();
      double sf = 2.0/m_gr->getL();
      m_dsp[N-1] = 0.0;
      m_dsp[N-2] =  2.0*(N-1)*sf*a[N-1]; 
      for(int k = N-2; k > 0; k--)
          m_dsp[k-1] = m_dsp[k+1] + 2.0*sf*k*a[k];
      std::copy(std::begin(m_dsp),std::end(m_dsp),std::begin(a));
  }

  void ChebX::dSX(const std::vector<double>& a,std::vector<double>& b,int n) 
  {
      std::copy(std::cbegin(a),std::cend(a),std::begin(b));
      dSX(b,n);
  }

  void ChebX::dSX(const std::vector<double>& a,std::vector<double>& b) 
  {
      std::copy(std::cbegin(a),std::cend(a),std::begin(b));
      dSX(b);
  }

}



//  void ChebX::dSX(std::vector<double>& a,int n)
//  {
//      if(n == 0)
//          return;
//      if(n < 0){ 
//          throw std::invalid_argument("ChebX::dSX(std::vector<double>& a,int n) \
//                  n must be a non-negative integer");
//      }
//      if(m_gr.gridType() == ChebGridType::EXTREMA){
//          int N = m_gr.getNx()-1;
//          double sf = 2.0/m_gr.getL();
//          //double sf = 1.0;
//          for(int i = 0; i < n; i++){
//              m_dico[N] = 0.0;
//              m_dico[N-1] =  2.0*N*a[N]*sf; // Weideman
//              for(int k = N-2; k > 0; k--)
//                m_dico[k] = (m_dico[k+2] + 2.0*sf*(k + 1)*a[k+1]) ;
//              m_dico[0] = 0.5*m_dico[2] + sf*a[1];
//              std::copy(std::begin(m_dico),std::end(m_dico),std::begin(a));
//          }
//      } else if(m_gr.gridType() == ChebGridType::ROOTS){
//          int N = m_gr.getNx()-1;
//          int nx = m_gr.getNx();
//          double sf = 2.0/m_gr.getL();
//          assert(m_gr.getL() == 2.0);
//          for(int i = 0; i < n; i++){
//              m_dico[N] = 0.0;
//              m_dico[N-1] =  2.0*N*a[N]*sf; // Weideman
//              for(int k = N-2; k > 0; k--)
//                m_dico[k] = (m_dico[k+2] + 2.0*sf*(k + 1)*a[k+1]) ;
//              m_dico[0] = 0.5*m_dico[2] + sf*a[1];
//              std::copy(std::begin(m_dico),std::end(m_dico),std::begin(a));
//          }
//
//          //for(int i = 0; i < n; i++){
//          //    m_dico[nx-1] = 0.0;
//          //    m_dico[nx-2] =  2.0*(nx-2)*a[nx-1]*sf; // Weideman
//          //    for(int k = nx-1; k > 1; k--)
//          //        m_dico[k-1] = m_dico[k+1] + 2.0*sf*(k - 1)*a[k];
//          //    std::copy(std::begin(m_dico),std::end(m_dico),std::begin(a));
//          //}
//      }
//  }

