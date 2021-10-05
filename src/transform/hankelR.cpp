
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <thread>
#include <boost/math/special_functions/bessel.hpp>
#include "spida/transform/hankelR.h"
#include "spida/grid/besselR.h" 

namespace spida {

  HankelTransformR::HankelTransformR(const BesselRootGridR& grid) : 
      m_nr(grid.getNr()),
      m_Ymk(grid.getNr()*grid.getNr())
  {
      m_alpha = grid.getjN()/pow(grid.getMaxSR(),2);
      initDHT(grid);
  }

  void HankelTransformR::R_To_SR(const double* in,double* out) 
  {
      for(unsigned m = 0; m < m_nr; m++){
          double sum = 0.0;
          for(unsigned k = 0; k < m_nr; k++)
              sum += m_Ymk[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
  }

  void HankelTransformR::R_To_SR(const dcmplx* in,dcmplx* out) 
  {
      for(unsigned m = 0; m < m_nr; m++){
          dcmplx sum = 0.0;
          for(unsigned k = 0; k < m_nr; k++)
              sum += m_Ymk[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
  }

  void HankelTransformR::SR_To_R(const double* in,double* out) 
  {
      for(unsigned k = 0; k < m_nr; k++){
          double sum = 0.0;
          for(unsigned m = 0; m < m_nr; m++)
              sum += m_Ymk[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
  }


  void HankelTransformR::SR_To_R(const dcmplx* in,dcmplx* out) 
  {
      for(unsigned k = 0; k < m_nr; k++){
          dcmplx sum = 0.0;
          for(unsigned m = 0; m < m_nr; m++)
              sum += m_Ymk[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
  }

  void HankelTransformR::initDHT(const BesselRootGridR& grid){
      const std::vector<double>& J0 = grid.getBesselRoots();
      std::vector<double> J1(m_nr);

      for(auto i = 0; i < m_nr; i++)
          J1[i] = boost::math::cyl_bessel_j<double>(1.0,J0[i]);

      double jN = grid.getjN();
      for(auto m = 0; m < m_nr; m++){
          for(auto k = 0; k < m_nr; k++){
              double beta_mk = 2.0/(jN*pow(J1[k],2));
              double arg = J0[m]*J0[k]/jN;
              double J0_mk = boost::math::cyl_bessel_j<double>(0.0,arg);
              m_Ymk[m*m_nr+k] = beta_mk*J0_mk;
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

      for(auto i = 0; i < m_nr; i++)
          J1[i] = boost::math::cyl_bessel_j<double>(1.0,J0[i]);

      double jN = grid.getjN();
      for(auto m = 0; m < m_nr; m++){
          for(auto k = 0; k < m_nr; k++){
              double beta_mk = 2.0/(jN*pow(J1[k],2));
              double arg = J0[m]*J0[k]/jN;
              double J0_mk = boost::math::cyl_bessel_j<double>(0.0,arg);
              m_Ymk[m*m_nr+k] = beta_mk*J0_mk;
          }
      }
  }






}










