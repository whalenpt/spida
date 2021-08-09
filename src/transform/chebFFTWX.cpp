
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <thread>
#include "spida/transform/chebFFTWX.h"
#include "spida/grid/chebX.h" 

namespace spida {

  ChebTransformFFTWX::ChebTransformFFTWX(const ChebExtremaGridX& grid) : 
      TransformX(grid),
      m_plan(new PlanChebTransformFFTWXdct(grid.getNx()))
  {
  }

  ChebTransformFFTWX::~ChebTransformFFTWX()
  {
      delete m_plan;
  }

  void ChebTransformFFTWX::X_To_SX(const std::vector<double>& in,std::vector<double>& out) 
  {
      m_plan->X_To_SX(in.data(),out.data());
  }

  void ChebTransformFFTWX::SX_To_X(const std::vector<double>& in,std::vector<double>& out)
  {
      m_plan->SX_To_X(in.data(),out.data());
  }

  // ChebTransformFFTWXdct flips [-1,1] x-grid to [0,pi] theta grid 
  PlanChebTransformFFTWXdct::PlanChebTransformFFTWXdct(int nx)
  {
    sz = nx;
    nlogic = 2*(nx-1);
    N = nx - 1;
    uFFT = new double[nx];
    for(int i = 0; i < nx; i++) uFFT[i] = 0.0;
    re_fftw_plan = fftw_plan_r2r_1d(nx,uFFT,uFFT,FFTW_REDFT00,FFTW_ESTIMATE);
    sp_fftw_plan = fftw_plan_r2r_1d(nx,uFFT,uFFT,FFTW_REDFT00,FFTW_ESTIMATE);
  }

  PlanChebTransformFFTWXdct::~PlanChebTransformFFTWXdct()
  {
    fftw_destroy_plan(re_fftw_plan);
    fftw_destroy_plan(sp_fftw_plan);
    delete [] uFFT;
  }

  void PlanChebTransformFFTWXdct::X_To_SX(const double* in,double* out) 
  {
    // flip [-1,1] x-grid to [0,pi] theta-grid
    for(int j = 0; j <= N; j++)  
      uFFT[N-j] = in[j]; 
    fftw_execute(re_fftw_plan);
    out[0] = uFFT[0]/(2.0*N);
    for (int j = 1; j < N; j++)  
      out[j] = uFFT[j]/N;
    out[N] = uFFT[N]/(2.0*N);
  }

  void PlanChebTransformFFTWXdct::SX_To_X(const double* in,double* out)  
  {
    uFFT[0] = in[0];
    for(int j = 1; j < N; j++)  
      uFFT[j] = in[j]/2.0; 
    uFFT[N] = in[N];
    fftw_execute(sp_fftw_plan);
    // flip [0,pi] theta-grid to [-1,1] x-grid
    for (int j = 0; j <= N; j++)  
      out[j] = uFFT[N-j];
  }

  void PlanChebTransformFFTWXdct::SpT_To_Re_2D(const double* in,double* out,int nd2)  
  {
    // flip [0,pi] theta-grid to [-1,1] x-grid
    uFFT[0] = in[0]; 
    for(int j = 1; j < N; j++) 
      uFFT[j] = in[j*nd2]/2.0; 
    uFFT[N] = in[N*nd2];
    fftw_execute(sp_fftw_plan);
    for (int j = 0; j < sz; j++)  
      out[j] = uFFT[sz-1-j];
  }

  void PlanChebTransformFFTWXdct::Re_To_SpT_2D(const double* in,double* out,int nd2) 
  {
    // flip [-1,1] x-grid to [0,pi] theta-grid
    for(int j = 0; j < sz; j++)  
      uFFT[sz-1-j] = in[j]; 
    fftw_execute(re_fftw_plan);

    out[0] = uFFT[0]/(2.0*N);
    for (int j = 1; j < N; j++)  
      out[j*nd2] = uFFT[j]/N;
    out[N*nd2] = uFFT[N]/(2.0*N);
  }

  void PlanChebTransformFFTWXdct::ReT_To_SpT_2D(const double* in,double* out,int nd2) 
  {
    // flip [-1,1] x-grid to [0,pi] theta-grid
    for(int j = 0; j < sz; j++)  
      uFFT[sz-1-j] = in[j*nd2]; 
    fftw_execute(re_fftw_plan);

    out[0] = uFFT[0]/(2.0*N);
    for (int j = 1; j < N; j++)  
      out[j*nd2] = uFFT[j]/N;
    out[N*nd2] = uFFT[N]/(2.0*N);
  }

  void PlanChebTransformFFTWXdct::SpT_To_ReT_2D(const double* in,double* out,int nd2)  
  {
    // flip [0,pi] theta-grid to [-1,1] x-grid
    uFFT[0] = in[0]; 
    for(int j = 1; j < N; j++) 
      uFFT[j] = in[j*nd2]/2.0; 
    uFFT[N] = in[N*nd2];
    fftw_execute(sp_fftw_plan);
    for (int j = 0; j < sz; j++)  
      out[j*nd2] = uFFT[sz-1-j];
  }


  PlanXr::PlanXr(int nx)
  {
    sz = nx;
    uFFT = new double[sz+2];
    for(int i = 0; i < sz+2; i++) uFFT[i] = 0.0;
    re_fftw_plan = fftw_plan_dft_r2c_1d(sz,uFFT,(fftw_complex*) uFFT, FFTW_ESTIMATE);
    sp_fftw_plan = fftw_plan_dft_c2r_1d(sz,(fftw_complex*) uFFT, uFFT, FFTW_ESTIMATE);
  }

  PlanXr::~PlanXr()
  {
    fftw_destroy_plan(re_fftw_plan);
    fftw_destroy_plan(sp_fftw_plan);
    delete [] uFFT;
  }

  void PlanXr::X_To_SX(const double* in,double* out)
  {
    for (int j = 0; j < sz; j++) 
      uFFT[j] = in[j]; 
    uFFT[sz] = 0.0;
    uFFT[sz+1] = 0.0;
    fftw_execute(re_fftw_plan);
    for (int j = 0; j <= sz/2; j++){  
      out[2*j] = uFFT[2*j];
      out[2*j+1] = uFFT[2*j+1];
    }
  }

  void PlanXr::SX_To_X(const double* in,double* out)
  {
    for(int j = 0; j < sz+2; j++) 
      uFFT[j] = in[j]/sz; 
    fftw_execute(sp_fftw_plan);
    for (int j = 0; j < sz; j++)  
      out[j] = uFFT[j];
  }


}

