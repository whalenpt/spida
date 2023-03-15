
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <thread>
#include <boost/math/special_functions/bessel.hpp>
#include "spida/transform/hankelR.h"
#include "spida/grid/besselR.h" 
#if defined(HAVE_OPENBLAS)
    #include "cblas.h"
#endif

namespace spida {

  HankelTransformR::HankelTransformR(const BesselRootGridR& grid) : 
      m_nr(grid.getNr()),
      m_Ymk(m_nr*m_nr),
      m_YmkC(m_nr*m_nr)
  {
      m_alpha = grid.getjN()/pow(grid.getMaxSR(),2);
      initDHT(grid);
  }

  void HankelTransformR::R_To_SR(const double* in,double* out) 
  {
/*
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE    {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
cblas_dgemv(const enum CBLAS_ORDER Order,
           const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
           const double alpha, const double *A, const int lda,
           const double *X, const int incX, const double beta,
           double *Y, const int incY);
dgemv y = alpha*A*x + beta*y
lda -> first dimension of A
*/

      #if defined(HAVE_OPENBLAS)
      cblas_dgemv(CblasRowMajor,CblasNoTrans,m_nr,m_nr,m_alpha,m_Ymk.data(),m_nr,in,1,0.0,out,1);
      #else
      for(unsigned m = 0; m < m_nr; m++){
          double sum = 0.0;
          for(unsigned k = 0; k < m_nr; k++)
              sum += m_Ymk[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
      #endif
  }

  void HankelTransformR::R_To_SR(const dcmplx* in,dcmplx* out) 
  {
      #if defined(HAVE_OPENBLAS)
      const double beta = 0.0;
      cblas_zgemv(CblasRowMajor,CblasNoTrans,m_nr,m_nr,&m_alpha,m_YmkC.data(),m_nr,in,1,&beta,out,1);
      #else
      for(unsigned m = 0; m < m_nr; m++){
          dcmplx sum = 0.0;
          for(unsigned k = 0; k < m_nr; k++)
              sum += m_YmkC[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
      #endif
  }

  void HankelTransformR::SR_To_R(const double* in,double* out) 
  {
      #if defined(HAVE_OPENBLAS)
      cblas_dgemv(CblasRowMajor,CblasNoTrans,m_nr,m_nr,1.0/m_alpha,m_Ymk.data(),m_nr,in,1,0.0,out,1);
      #else
      for(unsigned k = 0; k < m_nr; k++){
          double sum = 0.0;
          for(unsigned m = 0; m < m_nr; m++)
              sum += m_Ymk[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
      #endif
  }

  void HankelTransformR::SR_To_R(const dcmplx* in,dcmplx* out) 
  {
      #if defined(HAVE_OPENBLAS)
      const dcmplx a = 1.0/m_alpha;
      const dcmplx beta = 0.0;
      cblas_zgemv(CblasRowMajor,CblasNoTrans,m_nr,m_nr,&a,m_YmkC.data(),m_nr,in,1,&beta,out,1);
      #else
      for(unsigned k = 0; k < m_nr; k++){
          dcmplx sum = 0.0;
          for(unsigned m = 0; m < m_nr; m++)
              sum += m_YmkC[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
      #endif
  }

  void HankelTransformR::initDHT(const BesselRootGridR& grid){
      const std::vector<double>& J0 = grid.getBesselRoots();
      std::vector<double> J1(m_nr);

      for(unsigned i = 0; i < m_nr; i++)
          J1[i] = boost::math::cyl_bessel_j<double>(1.0,J0[i]);

      double jN = grid.getjN();
      for(unsigned m = 0; m < m_nr; m++){
          for(unsigned k = 0; k < m_nr; k++){
              double beta_mk = 2.0/(jN*pow(J1[k],2));
              double arg = J0[m]*J0[k]/jN;
              double J0_mk = boost::math::cyl_bessel_j<double>(0.0,arg);
              m_Ymk[m*m_nr+k] = beta_mk*J0_mk;
              m_YmkC[m*m_nr+k] = beta_mk*J0_mk;
          }
      }
  }


  HankelTransformRb::HankelTransformRb(const BesselRootGridR& grid,unsigned threads) : 
      m_threads(threads),
      m_nr(grid.getNr()),
      m_Ymk(grid.getNr()*grid.getNr())
  {
      m_alpha = grid.getjN()/pow(grid.getMaxSR(),2);
      initDHT(grid);
  }

  void HankelTransformRb::R_To_SR(const double* in,double* out) 
  {
      for(unsigned m = 0; m < m_nr; m++){
          double sum = 0.0;
          for(unsigned k = 0; k < m_nr; k++)
              sum += m_Ymk[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
  }

  void HankelTransformRb::R_To_SR(const dcmplx* in,dcmplx* out) 
  {
      std::vector<std::thread> workers;
      for(unsigned tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread([](\
                          unsigned tid,\
                          unsigned nthreads,\
                          unsigned nr,\
                          std::vector<double>& Ymk,\
                          double alpha,\
                          const dcmplx* v,\
                          dcmplx* w){
              for(unsigned m = tid; m < nr; m+=nthreads){
                  dcmplx sum = 0.0;
                  for(unsigned k = 0; k < nr; k++)
                      sum += Ymk[m*nr+k]*v[k];
                  w[m] = alpha*sum;
              }
          },tid,m_threads,m_nr,std::ref(m_Ymk),m_alpha,in,out));
      }

      for(auto& worker : workers){
          worker.join();
      }
  }

  void HankelTransformRb::SR_To_R(const double* in,double* out) 
  {
      for(unsigned k = 0; k < m_nr; k++){
          double sum = 0.0;
          for(unsigned m = 0; m < m_nr; m++)
              sum += m_Ymk[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
  }


  void HankelTransformRb::SR_To_R(const dcmplx* in,dcmplx* out) 
  {
      std::vector<std::thread> workers;
      for(unsigned tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread([](\
                          unsigned tid,\
                          unsigned nthreads,\
                          unsigned nr,\
                          std::vector<double>& Ymk,\
                          double alpha,\
                          const dcmplx* v,\
                          dcmplx* w){
              for(unsigned k = tid; k < nr; k+=nthreads){
                  dcmplx sum = 0.0;
                  for(unsigned m = 0; m < nr; m++)
                      sum += Ymk[k*nr+m]*v[m];
                  w[k] = sum/alpha;
              }
          },tid,m_threads,m_nr,std::ref(m_Ymk),m_alpha,in,out));
      }

      for(auto& worker : workers){
          worker.join();
      }
  }

  void HankelTransformRb::initDHT(const BesselRootGridR& grid){
      const std::vector<double>& J0 = grid.getBesselRoots();
      std::vector<double> J1(m_nr);

      for(unsigned i = 0; i < m_nr; i++)
          J1[i] = boost::math::cyl_bessel_j<double>(1.0,J0[i]);

      double jN = grid.getjN();
      for(unsigned m = 0; m < m_nr; m++){
          for(unsigned k = 0; k < m_nr; k++){
              double beta_mk = 2.0/(jN*pow(J1[k],2));
              double arg = J0[m]*J0[k]/jN;
              double J0_mk = boost::math::cyl_bessel_j<double>(0.0,arg);
              m_Ymk[m*m_nr+k] = beta_mk*J0_mk;
          }
      }
  }

}