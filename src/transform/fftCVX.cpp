#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "kiss_fft.h"
#include "spida/grid/uniformCVX.h" 
#include "spida/helper/constants.h"
#include "spida/transform/fftCVX.h"

namespace spida{

// kiss fft defined for x in [0,1], kx in [0,1]    
// FFTCVX wrapper will compute fft for x in [-xmin,xmax]  and kx in [-pi/L,pi/L] where L = xmax-xmin
// FFT{f(x)} = L*exp(-i*xmin*kx)*kissfft{f(x)}
FFTCVX::FFTCVX(const UniformGridCVX& grid) :
    m_nx(grid.getNx()),
    m_minx(grid.getMinX()),
    m_L(grid.getLX()),
    m_kx(grid.getSX()),
    m_temp(grid.getNsx())
{
    if(!((m_nx%2)==0))
        throw std::invalid_argument("Kiss fft requires even integer size");
    m_cfg_forward = kiss_fft_alloc(m_nx,0,nullptr,nullptr); 
    m_cfg_reverse = kiss_fft_alloc(m_nx,1,nullptr,nullptr);
}

FFTCVX::~FFTCVX()
{
    kiss_fft_free(m_cfg_forward);
    kiss_fft_free(m_cfg_reverse);
}

void FFTCVX::X_To_SX(const dcmplx* in,dcmplx* out) noexcept
{
    // Stuck with kiss_fft_cpx definition of complex -> use reinterpret and assume 
    // kiss interface won't change :(
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
    // Divide by FFT multiplier m_nx
    // Adjust phase since physical grid is not assumed to start at 0
    // Additional adjustment for length of grid (L)
    for(auto i = 0; i < m_nx; i++)
        out[i] *= (m_L*exp(ii*m_kx[i]*m_minx))/static_cast<double>(m_nx);
}
void FFTCVX::SX_To_X(const dcmplx* in,dcmplx* out) noexcept
{
    // Undo phase adjustment and length adjustment for inverse
    for(auto i = 0; i < m_nx; i++)
        m_temp[i] = in[i]*exp(-ii*m_kx[i]*m_minx)/m_L;

    // Stuck with kiss_fft_cpx definition of complex -> use reinterpret and assume 
    // kiss interface won't change :(
    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_temp.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
}

void FFTCVX::X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept { 
    X_To_SX(in.data(),out.data()); 
}

void FFTCVX::SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept { 
    SX_To_X(in.data(),out.data()); 
}
    
}