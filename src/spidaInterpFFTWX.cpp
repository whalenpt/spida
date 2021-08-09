
#include "spida/constants.h"
#include "spida/spidaInterpFFTWX.h"
#include "spida/spidaFFTWX.h"
#include "spida/grid/chebX.h"
#include "spida/interp.h"
#include <vector>
#include <memory>
#include <iostream>

namespace spida{

  ChebInterpFFTWX::ChebInterpFFTWX(const ChebExtremaGridX& grid) :
      m_gr(new ChebExtremaGridX(grid)),
      m_chebfftwx(new ChebFFTWX(grid)),
      m_ycheb(grid.getNx()),
      m_dycheb(grid.getNx())
  {
  }

  ChebInterpFFTWX::~ChebInterpFFTWX()
  {
      delete m_gr;
      delete m_chebfftwx;
  }


  void ChebInterpFFTWX::interpVec(const std::vector<double>& x_cheb,\
          const std::vector<double>& y_cheb,const std::vector<double>& x,\
          std::vector<double>& y)
  {
      SplineInterp spline(x_cheb,y_cheb);
      spline.eval(x,y);
  }

  double ChebInterpFFTWX::interpVal(const std::vector<double>& x_cheb,\
          const std::vector<double>& y_cheb,double xout)
  {
      SplineInterp spline(x_cheb,y_cheb);
      return spline.eval(xout);
  }

  double ChebInterpFFTWX::dXInterp(const std::vector<double>& x,\
          const std::vector<double>& y,double xout,int n) 
  {
      auto minIndx = 0;
      auto maxIndx = 0;
      auto refIndx = 0;
      for(auto i = 0; i < x.size(); i++){
        if(x[i] < xout)
          refIndx++;
        else
          break;
      }

      int nx = m_gr->getNx();
      if(refIndx + nx/2 < x.size())
        maxIndx = refIndx + nx/2;
      else
        maxIndx = x.size() - 1; 
      if(refIndx - nx/2 >= 0)
        minIndx = refIndx - nx/2;
      else
        minIndx = 0; 

      interpVec(x,y,m_gr->getX(),m_ycheb);
      m_chebfftwx->dX(m_ycheb,m_dycheb,n);
      return interpVal(m_gr->getX(),m_dycheb,xout);
  }

  void ChebInterpFFTWX::dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
          const std::vector<double>& xout,std::vector<double>& dyout,int n)
  {
      interpVec(xin,yin,m_gr->getX(),m_ycheb);
      m_chebfftwx->dX(m_ycheb,m_dycheb,n);
      interpVec(m_gr->getX(),m_dycheb,xout,dyout);
  }
    
  void ChebInterpFFTWX::dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
          std::vector<double>& dyout,int n) 
  {
      interpVec(xin,yin,m_gr->getX(),m_ycheb);
      m_chebfftwx->dX(m_ycheb,m_dycheb,n);
      interpVec(m_gr->getX(),m_dycheb,xin,dyout);
  }

//  void ChebInterpFFTWX::iX(const std::vector<double>& x,const std::vector<double>& y,std::vector<double>& yout,int n) 
//  {
//    spida::ChebGridX gr(m_nx,x[0],x.back());
//    m_L = gr.getL();
//    interpVec(x,y,gr.getX(),m_ycheb);
//    m_tr->X_To_SX(m_ycheb,m_sp);
//    iSX(m_sp,n);
//    m_tr->SX_To_X(m_sp,m_ycheb);
//    interpVec(gr.getX(),m_ycheb,x,yout);
//  }
//
//  double ChebInterpFFTWX::iX(const std::vector<double>& x,const std::vector<double>& y,double xmin,double xmax,int n) 
//  {
//    spida::ChebGridX gr(m_nx,x[0],x.back());
//    m_L = gr.getL();
//    interpVec(x,y,gr.getX(),m_ycheb);
//    m_tr->X_To_SX(m_ycheb,m_sp);
//    iSX(m_sp,n);
//    m_tr->SX_To_X(m_sp,m_dico);
//    int minI = gr.getIndexX(xmin);
//    int maxI = gr.getIndexX(xmax);
//    return (m_dico[maxI] - m_dico[minI]);
//  }

//  void ChebInterpFFTWX::iSX(const std::vector<double>& a,std::vector<double>& b,int n) 
//  {
//    if(n == 0){
//      for(int i = 0; i < m_nx; i++)
//        b[i] = a[i];
//    }
//    else if(n == 1){
//      iS1(a,b);
//    }
//    else if(n > 1){
//      iS1(a,b);
//      int diffCount = 1;
//      while(diffCount < n){
//        iS1(b);
//        diffCount++;
//      }
//    }
//    else{
//        throw pw::Exception("ChebInterpFFTWX::iSX(const std::vector<double>& a,std::vector<double>& b,int n)","n must be"\
//                " a non-negative integer");
//    }
//  }
//
//  void ChebInterpFFTWX::iSX(std::vector<double>& a,int n)
//  {
//    if(n == 1)
//      iS1(a);
//    else if(n > 1){
//      int diffCount = 0;
//      while(diffCount < n){
//        iS1(a);
//        diffCount++;
//      }
//    }
//    else if(n == 0)
//      return;
//    else{
//        throw pw::Exception("ChebInterpFFTWX::iSX(std::vector<double>& a,int n)","n must be"\
//                " a non-negative integer");
//    }
//  }
//
//  void ChebInterpFFTWX::iS1(std::vector<double>& a) 
//  {
//    int N = m_nx-1;
//    double sf = m_L/2.0;
//    m_dico[0] = a[1]/4.0;
//    m_dico[1] = a[0] - a[2]/2.0;
//    for(int k = 2; k < N; k++)
//      m_dico[k] = (a[k-1] - a[k+1])/(2*k);
//    m_dico[N] = a[N-1]/(2.0*N);
//    for(int i = 0; i <= N; i++)
//      a[i] = sf*m_dico[i];
//  }
//
//  void ChebInterpFFTWX::iS1(const std::vector<double>& a,std::vector<double>& b)  
//  {
//    int N = m_nx-1;
//    double sf = m_L/2.0; // y = xmin + L*(1+x)/2 -> dx = (2/L)dy
//    b[0] = sf*a[1]/4.0;
//    b[1] = sf*(a[0] - a[2]/2.0);
//    for(int k = 2; k < N; k++)
//      b[k] = sf*((a[k-1] - a[k+1])/(2*k));
//    b[N] = sf*(a[N-1]/(2.0*N));
//  }


}




