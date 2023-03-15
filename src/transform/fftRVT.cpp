#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "spida/grid/uniformRVT.h" 
#include "spida/helper/constants.h"
#include "spida/transform/fftRVT.h"

namespace spida{

FFTRVT::FFTRVT(const UniformGridRVT& grid) :
    m_nt(grid.getNt()),
    m_nst(grid.getNst()),
    m_minI(grid.getMinI()),
    m_maxI(grid.getMaxI()),
    m_rFFTr(grid.getNt(),0.0),
    m_temp(grid.getNt(),0.0),
    m_cFFT(grid.getNt(),0.0),
    m_omega(grid.getST()),
    m_mint(grid.getMinT()),
    m_L(grid.getLT())
{
    if(!((m_nt%2)==0))
        throw std::invalid_argument("Kiss real fft requires even integer size");

    // T-transform centered around -iwt -> use inverse kissfft for forward direction (s.t.
    // fft central frequency is the first half of the grid and aligns with the UniformGridT class)
    m_cfg_forward = kiss_fft_alloc(m_nt,1,nullptr,nullptr); 
    m_cfg_reverse = kiss_fft_alloc(m_nt,0,nullptr,nullptr);

    m_rcfg_forward = kiss_fftr_alloc(m_nt,0,nullptr,nullptr);   
    m_rcfg_reverse = kiss_fftr_alloc(m_nt,1,nullptr,nullptr);   
}

FFTRVT::~FFTRVT()
{
    kiss_fft_free(m_cfg_forward);
    kiss_fft_free(m_cfg_reverse);
    kiss_fft_free(m_rcfg_forward);
    kiss_fft_free(m_rcfg_reverse);

}

void FFTRVT::T_To_ST(const double* in,dcmplx* out)
{
    kiss_fftr(m_rcfg_forward,reinterpret_cast<const kiss_fft_scalar*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(m_temp.data()));
    // Divide by FFT multiplier m_nt and adjust phase since physical grid is not assumed to start at 0
    // Also take negative imaginary component -> (-iwt convention vs kiss fft convention of iwt)
    for(auto j = m_minI; j <= m_maxI; j++){
        dcmplx tempval((m_L*exp(ii*m_omega[j-m_minI]*m_mint))*m_temp[j]/static_cast<double>(m_nt));
        out[j-m_minI] = dcmplx(tempval.real(),-tempval.imag());
    }
}

void FFTRVT::ST_To_T(const dcmplx* in,double* out)
{
    // band limited minimum frequency 
    for(unsigned j = 0; j < m_minI; j++) 
        m_temp[j] = 0.0; 
    for(unsigned j = m_minI; j <= m_maxI; j++){
        dcmplx tempval(exp(-ii*m_omega[j-m_minI]*m_mint)*in[j-m_minI]/m_L);
        m_temp[j] = dcmplx(tempval.real(),-tempval.imag());
    }
    // band limited maximum frequency
    for(unsigned j = m_maxI+1; j < m_nt/2+1; j++)
        m_temp[j] = 0.0; 
    kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_temp.data()),\
                  reinterpret_cast<kiss_fft_scalar*>(out));
}

void FFTRVT::CVT_To_ST(const dcmplx* in,dcmplx* out)
{
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in),\
            reinterpret_cast<kiss_fft_cpx*>(m_cFFT.data()));
    for(unsigned j = m_minI; j<=m_maxI; j++)
        out[j-m_minI] = m_L*exp(ii*m_omega[j-m_minI]*m_mint)*m_cFFT[j]/static_cast<double>(m_nt);
}

void FFTRVT::ST_To_CVT(const dcmplx* in,dcmplx* out)
{
    std::fill(m_cFFT.begin(),m_cFFT.begin()+m_minI,dcmplx(0,0));
    for(auto j = m_minI; j <= m_maxI; j++) 
        m_cFFT[j] = exp(-ii*m_omega[j-m_minI]*m_mint)*in[j-m_minI]/m_L;
    std::fill(m_cFFT.begin()+m_maxI+1,m_cFFT.end(),dcmplx(0,0));
    kiss_fft(m_cfg_reverse,reinterpret_cast<kiss_fft_cpx*>(m_cFFT.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out));
}

void FFTRVT::T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out) { 
    T_To_ST(in.data(),out.data()); 
}

void FFTRVT::ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out) { 
    ST_To_T(in.data(),out.data()); 
}

void FFTRVT::CVT_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) { 
    CVT_To_ST(in.data(),out.data()); 
}

void FFTRVT::ST_To_CVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) { 
    ST_To_CVT(in.data(),out.data()); 
}


}