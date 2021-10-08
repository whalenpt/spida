#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include "spida/transform/fftRVX.h"
#include "spida/grid/uniformRVX.h" 
#include "spida/helper/constants.h"
#include "kiss_fftr.h"

namespace spida{

FFTRVX::FFTRVX(const UniformGridRVX& grid) :
    m_nx(grid.getNx()),
    m_minx(grid.getMinX()),
    m_L(grid.getLX()),
    m_kx(grid.getSX()),
    m_temp(grid.getNx()/2+1)
{
    if(!((m_nx%2)==0))
        throw std::invalid_argument("Kiss fftr requires even integer size");
    m_rcfg_forward = kiss_fftr_alloc(m_nx,0,nullptr,nullptr);   
    m_rcfg_reverse = kiss_fftr_alloc(m_nx,1,nullptr,nullptr);   
}

FFTRVX::~FFTRVX()
{
    kiss_fft_free(m_rcfg_forward);
    kiss_fft_free(m_rcfg_reverse);
}

void FFTRVX::X_To_SX(const double* in,dcmplx* out) noexcept
{
    kiss_fftr(m_rcfg_forward,reinterpret_cast<const kiss_fft_scalar*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
    // Divide by FFT multiplier m_nx and adjust phase since physical grid is not assumed to start at 0
    for(auto i = 0; i < m_nx/2+1; i++)
        out[i] *= (m_L*exp(ii*m_kx[i]*m_minx))/static_cast<double>(m_nx);
}

void FFTRVX::SX_To_X(const dcmplx* in,double* out) noexcept
{
    // Undo phase adjustment for inverse
    for(auto i = 0; i < m_nx/2+1; i++)
        m_temp[i] = in[i]*(exp(-ii*m_kx[i]*m_minx)/m_L);
    kiss_fftri(m_rcfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_temp.data()),\
                  reinterpret_cast<kiss_fft_scalar*>(out));
}

void FFTRVX::X_To_SX(const std::vector<double>& in,std::vector<dcmplx>& out) noexcept
{ X_To_SX(in.data(),out.data()); }


void FFTRVX::SX_To_X(const std::vector<dcmplx>& in,std::vector<double>& out) noexcept
{ SX_To_X(in.data(),out.data()); }



}







