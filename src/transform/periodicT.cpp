
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include "spida/transform/periodicT.h"
#include "spida/grid/uniformT.h" 
#include "spida/constants.h"
#include "kiss_fft.h"
#include "kiss_fftr.h"

namespace spida{

PeriodicTransformT::PeriodicTransformT(const UniformGridT& grid) :
    m_rFFTr(grid.getNt(),0.0),
    m_rFFTs(grid.getNt()/2+1,0.0),
    m_cFFT(grid.getNt(),0.0)
{
    m_nt = grid.getNt();
    if(!((m_nt%2)==0))
        throw std::invalid_argument("Kiss real fft requires even integer size");

    m_nst = grid.getNst();
    m_minI = grid.getMinI();
    m_maxI = grid.getMaxI();
    // T-transform centered around -iwt -> use inverse kissfft for forward direction (s.t.
    // fft central frequency is the first half of the grid and aligns with the UniformGridT class)
    m_cfg_forward = kiss_fft_alloc(m_nt,1,nullptr,nullptr); 
    m_cfg_reverse = kiss_fft_alloc(m_nt,0,nullptr,nullptr);

    // Can't use inverse kissfft for forward direction with real fields (different memory allocation sizes)
    m_rcfg_forward = kiss_fftr_alloc(m_nt,0,nullptr,nullptr);   
    m_rcfg_reverse = kiss_fftr_alloc(m_nt,1,nullptr,nullptr);   
}


PeriodicTransformT::~PeriodicTransformT()
{
    kiss_fft_free(m_cfg_forward);
    kiss_fft_free(m_cfg_reverse);
    kiss_fft_free(m_rcfg_forward);
    kiss_fft_free(m_rcfg_reverse);

}

void PeriodicTransformT::execute_forward()
{
    kiss_fftr(m_rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(m_rFFTr.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(m_rFFTs.data()));
}

void PeriodicTransformT::execute_backward()
{
    kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_rFFTs.data()),\
                  reinterpret_cast<kiss_fft_scalar*>(m_rFFTr.data()));
}

void PeriodicTransformT::T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out)
{
    // FFTW FORWARD - > +iwt transform def -> reverse time t -> -t
    std::reverse_copy(std::begin(in),std::end(in),std::begin(m_rFFTr));
    execute_forward();
    std::copy(std::begin(m_rFFTs)+m_minI,std::begin(m_rFFTs)+m_maxI+1,std::begin(out));
}

void PeriodicTransformT::T_To_ST(const double* in,dcmplx* out)
{
    // FFTW FORWARD - > +iwt transform def -> reverse time t -> -t
    for (unsigned int j = 0; j < m_nt; j++) 
        m_rFFTr[j] = in[m_nt-j-1]; 
    execute_forward();
    for(unsigned int j = m_minI; j<=m_maxI; j++)
        out[j-m_minI] = m_rFFTs[j];
}

void PeriodicTransformT::ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out)
{
    std::fill(m_rFFTs.begin(),m_rFFTs.begin()+m_minI,dcmplx(0,0));
    for(auto j = m_minI; j <= m_maxI; j++)
        m_rFFTs[j] = in[j-m_minI]/static_cast<double>(m_nt); 
    std::fill(m_rFFTs.begin()+m_maxI+1,m_rFFTs.end(),dcmplx(0,0));
    execute_backward();
    std::reverse_copy(std::begin(m_rFFTr),std::end(m_rFFTr),std::begin(out));
}

void PeriodicTransformT::ST_To_T(const dcmplx* in,double* out)
{
    // band limited minimum frequency 
    for(unsigned int j = 0; j < m_minI; j++) 
        m_rFFTs[j] = 0.0; 
    for(unsigned int j = m_minI; j <= m_maxI; j++)
        m_rFFTs[j] = in[j-m_minI]/static_cast<double>(m_nt); 
    // band limited maximum frequency
    for(unsigned int j = m_maxI+1; j < m_nt/2+1; j++)
        m_rFFTs[j] = 0.0; 
    execute_backward();
    // FFTW BACKWARD - > -iwt transform def -> reverse time t -> -t
    for (unsigned int j = 0; j < m_nt; j++)  
      out[j] = m_rFFTr[m_nt-j-1];
}


void PeriodicTransformT::T_To_ST_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(m_cFFT.data()));
    std::copy(std::begin(m_cFFT)+m_minI,std::begin(m_cFFT)+m_maxI+1,std::begin(out));
}

void PeriodicTransformT::ST_To_T_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    std::fill(m_cFFT.begin(),m_cFFT.begin()+m_minI,dcmplx(0,0));
    for(auto j = m_minI; j <= m_maxI; j++) 
        m_cFFT[j] = in[j-m_minI]/static_cast<double>(m_nt); 
    // Set negative half of spectrum to zero (fine for real fields)
    std::fill(m_cFFT.begin()+m_maxI+1,m_cFFT.end(),dcmplx(0,0));
    kiss_fft(m_cfg_reverse,reinterpret_cast<kiss_fft_cpx*>(m_cFFT.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
}



}


