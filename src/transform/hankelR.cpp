
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>
//#include <gsl/gsl_sf_bessel.h>
#include "spida/transform/hankelR.h"
#include "spida/grid/besselR.h" 

namespace spida {

  void printHankel(const HankelTransformR& transform,std::ostream& os)
  {
      const std::vector<double>& Ymk = transform.getYmk();
      int nr = transform.getNr();
      os << std::endl << std::endl;
      for(auto m = 0; m < nr; m++){
          for(auto k = 0; k < nr-1; k++)
              os << std::setprecision(3) << std::scientific << Ymk[m*nr+k] << " ";
          os << std::setprecision(3) << std::scientific << Ymk[m*nr+nr-1]  << std::endl;
      }
      os << std::endl;
  }

  HankelTransformR::HankelTransformR(const BesselRootGridR& grid) : 
      TransformR(grid), 
      m_nr(grid.getNr()),
      m_Ymk(grid.getNr()*grid.getNr())
  {
      m_alpha = grid.getjN()/pow(grid.getMaxSR(),2);
      initDHT(grid);
  }

  void HankelTransformR::R_To_SR(const double* in,double* out) 
  {
      for(unsigned int m = 0; m < m_nr; m++){
          double sum = 0.0;
          for(unsigned int k = 0; k < m_nr; k++)
              sum += m_Ymk[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
  }

  void HankelTransformR::R_To_SR(const dcmplx* in,dcmplx* out) 
  {
      for(unsigned int m = 0; m < m_nr; m++){
          dcmplx sum = 0.0;
          for(unsigned int k = 0; k < m_nr; k++)
              sum += m_Ymk[m*m_nr+k]*in[k];
          out[m] = m_alpha*sum;
      }
  }

  void HankelTransformR::SR_To_R(const double* in,double* out) 
  {
      for(unsigned int k = 0; k < m_nr; k++){
          double sum = 0.0;
          for(unsigned int m = 0; m < m_nr; m++)
              sum += m_Ymk[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
  }


  void HankelTransformR::SR_To_R(const dcmplx* in,dcmplx* out) 
  {
      for(unsigned int k = 0; k < m_nr; k++){
          dcmplx sum = 0.0;
          for(unsigned int m = 0; m < m_nr; m++)
              sum += m_Ymk[k*m_nr+m]*in[m];
          out[k] = sum/m_alpha;
      }
  }

  void HankelTransformR::initDHT(const BesselRootGridR& grid){
      const std::vector<double>& J0 = grid.getBesselRoots();
      std::vector<double> J1(m_nr);

      for(auto i = 0; i < m_nr; i++)
          J1[i] = boost::math::cyl_bessel_j<double>(1.0,J0[i]);

//      for(auto i = 0; i < m_nr; i++)
//          J1[i] = gsl_sf_bessel_J1(J0[i]);

      double jN = grid.getjN();
      for(auto m = 0; m < m_nr; m++){
          for(auto k = 0; k < m_nr; k++){
              double beta_mk = 2.0/(jN*pow(J1[k],2));
              double arg = J0[m]*J0[k]/jN;
              double J0_mk = boost::math::cyl_bessel_j<double>(0.0,arg);
              //double J0_mk = gsl_sf_bessel_J0(arg);
              m_Ymk[m*m_nr+k] = beta_mk*J0_mk;
          }
      }
  }
}










