#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include "spida/transform/fftX.h"
#include "spida/grid/uniformX.h" 
#include "spida/helper/constants.h"
#include "kiss_fft.h"

namespace spida{

FFTX::FFTX(const UniformGridX& grid) 
{
    m_nx = grid.getNx();
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

void FFTX::X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out.data()));
}

void FFTX::X_To_SX(const dcmplx* in,dcmplx* out)
{
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
}

void FFTX::SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(in.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out.data()));
    for(auto i = 0; i < m_nx; i++)
        out[i] /= static_cast<double>(m_nx);
}

void FFTX::SX_To_X(const dcmplx* in,dcmplx* out)
{
    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
    for(auto i = 0; i < m_nx; i++)
        out[i] /= static_cast<double>(m_nx);
}


}







