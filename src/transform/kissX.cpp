
#include <fftw3.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <thread>
#include "spida/transform/transformX.h"
#include "spida/grid/gridX.h" 

namespace spida {

  // ChebTransformXdct flips [-1,1] x-grid to [0,pi] theta grid 
  PlanChebTransformXdct::PlanChebTransformXdct(int nx) : m_uFFT(nx,0.0),
    m_uFFTs(nx,0.0)
  {
      m_nx = nx;
      m_N = nx - 1;
      m_nlogic = 2*(nx-1);
      m_rcfg_forward = kiss_fftr_alloc(m_nx,0,nullptr,nullptr);   
      m_rcfg_reverse = kiss_fftr_alloc(m_nx,1,nullptr,nullptr);   
  }

  PlanChebTransformXdct::~PlanChebTransformXdct()
  {
      kiss_fft_free(m_rcfg_forward);
      kiss_fft_free(m_rcfg_reverse);
  }

  void PlanChebTransformXdct::X_To_SX(const std::vector<double>& in,std::vector<double>& out) 
  {
    // flip [-1,1] x-grid to [0,pi] theta-grid
    std::reverse_copy(std::begin(in),std::end(in),std::begin(m_uFFT));
    kiss_fftr(m_rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(m_uFFTs.data()));
    for (int j = 0; j < m_N; j++)  
        out[j] = m_uFFTs[j].real()/m_N;
  }

  void PlanChebTransformXdct::SX_To_X(const std::vector<double>& in,std::vector<double>& out)  
  {
    // flip [0,pi] theta-grid to [-1,1] x-grid
    for (int j = 0; j < m_nx; j++)  
        m_uFFTs[j] = in[j];

    kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_uFFTs.data()),\
                  reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()));
    std::reverse_copy(std::begin(m_uFFT),std::end(m_uFFT),std::begin(out));
  }

  void PlanChebTransformXdct::SpT_To_Re_2D(const std::vector<double>& in,std::vector<double>& out,int nd2)  
  {
    // flip [0,pi] theta-grid to [-1,1] x-grid
    m_uFFTs[0] = in[0]; 
    for(int j = 0; j < m_nx; j++) 
      m_uFFTs[j] = in[j*nd2]/2.0; 
    uFFTs[m_N] = in[m_N*nd2];

    kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_uFFTs.data()),\
                  reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()));
    std::reverse_copy(std::begin(m_uFFT),std::end(m_uFFT),std::begin(out));
  }

  void PlanChebTransformXdct::Re_To_SpT_2D(const std::vector<double>& in,std::vector<double>& out,int nd2) 
  {
    // flip [-1,1] x-grid to [0,pi] theta-grid
    std::reverse_copy(std::begin(in),std::end(in),std::begin(uFFT));
    kiss_fftr(m_rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(m_uFFTs.data()));

    out[0] = m_uFFTs[0].real()/(2.0*m_N);
    for (int j = 1; j < m_N; j++)  
      out[j*nd2] = m_uFFTs[j].real()/m_N;
    out[N*nd2] = m_uFFTs[m_N].real()/(2.0*m_N);
  }

  void PlanChebTransformXdct::ReT_To_SpT_2D(const std::vector<double>& in,std::vector<double>& out,int nd2) 
  {
    // flip [-1,1] x-grid to [0,pi] theta-grid
    for(int j = 0; j < m_nx; j++)  
      m_uFFT[m_nx-1-j] = in[j*nd2]; 

    kiss_fftr(m_rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(m_uFFTs.data()));
    out[0] = uFFTs[0].real()/(2.0*m_N);
    for (int j = 1; j < N; j++)  
      out[j*nd2] = m_uFFTs[j].real()/m_N;
    out[m_N*nd2] = m_uFFTs[N].real()/(2.0*m_N);
  }

  void PlanChebTransformXdct::SpT_To_ReT_2D(const std::vector<double>& in,std::vector<double>& out,int nd2)  
  {
    // flip [0,pi] theta-grid to [-1,1] x-grid
    m_uFFTs[0] = in[0]; 
    for(int j = 1; j < m_N; j++) 
      m_uFFTs[j] = in[j*nd2]/2.0; 
    m_uFFTs[m_N] = in[m_N*nd2];
    kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_uFFTs.data()),\
                  reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()));
    for (int j = 0; j < m_nx; j++)  
      out[j*nd2] = m_uFFT[m_nx-1-j];
  }


  ChebTransformX::ChebTransformX(int cnx) 
  {
    m_plan = new PlanChebTransformXdct(cnx);
  }

  ChebTransformX::~ChebTransformX()
  {
    delete m_plan;
  }

  void ChebTransformX::X_To_SX(const std::vector<double>& in,std::vector<double>& out) 
  {
    m_plan->X_To_SX(&in[0],&out[0]);
  }

  void ChebTransformX::SX_To_X(const std::vector<double>& in,std::vector<double>& out)
  {
    m_plan->SX_To_X(&in[0],&out[0]);
  }

  PlanXr::PlanXr(int nx) : m_nx(nx),
    m_uFFT(nx,0.0), 
    m_uFFTs(nx,0.0)
  {
      m_rcfg_forward = kiss_fftr_alloc(m_nx,0,nullptr,nullptr);   
      m_rcfg_reverse = kiss_fftr_alloc(m_nx,1,nullptr,nullptr);   
  }

  PlanXr::~PlanXr()
  {
      kiss_fft_free(m_rcfg_forward);
      kiss_fft_free(m_rcfg_reverse);
  }

  void PlanXr::X_To_SX(const std::vector<double>& in,std::vector<dcmplx>& out)
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_uFFT));
      kiss_fftr(m_rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(m_uFFTs.data()));
      std::copy(std::begin(m_uFFTs),std::end(m_uFFTs),std::begin(out));
  }

  void PlanXr::SX_To_X(const std::vector<dcmplx>& in,std::vector<double>& out)
  {
      for(auto j = 0; j < m_nx; j++)
          m_rFFTs[j] = in[j]/static_cast<double>(m_nx); 
      kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_uFFTs.data()),\
                    reinterpret_cast<kiss_fft_scalar*>(m_uFFT.data()));
      std::copy(std::begin(m_uFFT),std::end(m_uFFT),std::begin(out));
  }
}

