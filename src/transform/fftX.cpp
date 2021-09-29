#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include "spida/transform/fftX.h"
#include "spida/grid/uniformX.h" 
#include "spida/helper/constants.h"
#include "kiss_fft.h"

namespace spida{

FFTX::FFTX(const UniformGridX& grid) :
    m_nx(grid.getNx()),
    m_minx(grid.getMinX()),
    m_L(grid.getLX()),
    m_kx(grid.getSX()),
    m_temp(grid.getNx())
{
    if(!((m_nx%2)==0))
        throw std::invalid_argument("Kiss fft requires even integer size");
    m_cfg_forward = kiss_fft_alloc(m_nx,0,nullptr,nullptr); 
    m_cfg_reverse = kiss_fft_alloc(m_nx,1,nullptr,nullptr);
}

FFTX::~FFTX()
{
    kiss_fft_free(m_cfg_forward);
    kiss_fft_free(m_cfg_reverse);
}

void FFTX::X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept
{
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out.data()));
    // Divide by FFT multiplier m_nx and adjust phase since physical grid is not assumed to start at 0
    for(auto i = 0; i < out.size(); i++)
        out[i] *= exp(ii*m_kx[i]*m_minx)/static_cast<double>(m_nx);
//    for(auto i = 0; i < out.size(); i++)
//        out[i] /= static_cast<double>(m_nx);

}

void FFTX::X_To_SX(const dcmplx* in,dcmplx* out) noexcept
{
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
    // Divide by FFT multiplier m_nx and adjust phase since physical grid is not assumed to start at 0
    for(auto i = 0; i < m_nx; i++)
        out[i] *= exp(ii*m_kx[i]*m_minx)/static_cast<double>(m_nx);

//    for(auto i = 0; i < m_nx; i++)
//        out[i] /= static_cast<double>(m_nx);
}

void FFTX::SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept
{
    // Undo phase adjustment for inverse
    for(auto i = 0; i < in.size(); i++)
        m_temp[i] = in[i]*exp(-ii*m_kx[i]*m_minx);
    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_temp.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out.data()));

//    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(in.data()),\
//                  reinterpret_cast<kiss_fft_cpx*>(out.data()));

}

void FFTX::SX_To_X(const dcmplx* in,dcmplx* out) noexcept
{
    // Undo phase adjustment for inverse
    for(auto i = 0; i < m_nx; i++)
        m_temp[i] = in[i]*exp(-ii*m_kx[i]*m_minx);

    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_temp.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out));

//    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(in),\
//                  reinterpret_cast<kiss_fft_cpx*>(out));

}



}







