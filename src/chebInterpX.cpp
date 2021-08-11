
#include "spida/constants.h"
#include "spida/spidaInterpX.h"
#include "spida/spidaChebX.h"
#include "spida/grid/chebX.h"
#include "spida/interp.h"
#include <vector>
#include <memory>
#include <iostream>

namespace spida{

  ChebInterpX::ChebInterpX(int ninterp,int minx,int maxx) :
      m_chebx(new ChebX(ChebRootGridX(ninterp,minx,maxx))),
      m_ycheb(ninterp),
      m_dycheb(ninterp)
  {
  }

  ChebInterpX::~ChebInterpX()
  {
      delete m_chebx;
  }

  double ChebInterpX::dXInterp(const std::vector<double>& xin,\
          const std::vector<double>& yin,double xout,int n) 
  {
      SplineInterp spline_to_cheb(xin,yin);
      spline_to_cheb.eval(m_chebx->getGridX().getX(),m_ycheb);
      m_chebx->dX(m_ycheb,m_dycheb,n);
      SplineInterp spline_from_cheb(m_chebx->getGridX().getX(),m_dycheb);
      return spline_from_cheb.eval(xout);
  }

  void ChebInterpX::dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
          const std::vector<double>& xout,std::vector<double>& dyout,int n)
  {
      SplineInterp spline_to_cheb(xin,yin);
      spline_to_cheb.eval(m_chebx->getGridX().getX(),m_ycheb);
      m_chebx->dX(m_ycheb,m_dycheb,n);
      SplineInterp spline_from_cheb(m_chebx->getGridX().getX(),m_dycheb);
      return spline_from_cheb.eval(xout,dyout);
  }
    
//  void ChebInterpX::iX(const std::vector<double>& x,const std::vector<double>& y,std::vector<double>& yout,int n) 
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
//  double ChebInterpX::iX(const std::vector<double>& x,const std::vector<double>& y,double xmin,double xmax,int n) 
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

//  void ChebInterpX::iSX(const std::vector<double>& a,std::vector<double>& b,int n) 
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
//        throw pw::Exception("ChebInterpX::iSX(const std::vector<double>& a,std::vector<double>& b,int n)","n must be"\
//                " a non-negative integer");
//    }
//  }
//
//  void ChebInterpX::iSX(std::vector<double>& a,int n)
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
//        throw pw::Exception("ChebInterpX::iSX(std::vector<double>& a,int n)","n must be"\
//                " a non-negative integer");
//    }
//  }
//
//  void ChebInterpX::iS1(std::vector<double>& a) 
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
//  void ChebInterpX::iS1(const std::vector<double>& a,std::vector<double>& b)  
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








