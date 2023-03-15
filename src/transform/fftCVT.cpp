#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "kiss_fft.h"
#include "spida/grid/uniformCVT.h" 
#include "spida/helper/constants.h"
#include "spida/transform/fftCVT.h"

namespace spida{

// kiss fft defined for t in [0,1], omega in [0,1]    
// FFTCVT wrapper will compute fft for t in [-tmin,tmax]  and omega in [-pi/L,pi/L] where L = tmax-tmin
// FFT{f(t)} = L*exp(i*tmin*omega)*kissfft{f(t)}
FFTCVT::FFTCVT(const UniformGridCVT& grid) :
    m_nt(grid.getNt()),
    m_mint(grid.getMinT()),
    m_L(grid.getLT()),
    m_omega(grid.getST()),
    m_temp(grid.getNst())
{
    if(!((m_nt%2)==0))
        throw std::invalid_argument("Kiss fft requires even integer size");

    // T-transform centered around -iwt -> use inverse kissfft for forward direction
    m_cfg_forward = kiss_fft_alloc(m_nt,1,nullptr,nullptr); 
    m_cfg_reverse = kiss_fft_alloc(m_nt,0,nullptr,nullptr);
}

FFTCVT::~FFTCVT()
{
    kiss_fft_free(m_cfg_forward);
    kiss_fft_free(m_cfg_reverse);
}

void FFTCVT::T_To_ST(const dcmplx* in,dcmplx* out) noexcept
{
    // Stuck with kiss_fft_cpx definition of complex -> use reinterpret and assume 
    // kiss interface won't change :(
    kiss_fft(m_cfg_forward,reinterpret_cast<const kiss_fft_cpx*>(in),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
    // Divide by FFT multiplier m_nx
    // Adjust phase since physical grid is not assumed to start at 0
    // Additional adjustment for length of grid (L)
    for(unsigned i = 0; i < m_nt; i++)
        out[i] *= (m_L*exp(-ii*m_omega[i]*m_mint))/static_cast<double>(m_nt);
}

void FFTCVT::ST_To_T(const dcmplx* in,dcmplx* out) noexcept
{
    // Undo phase adjustment and length adjustment for inverse
    for(unsigned i = 0; i < m_nt; i++)
        m_temp[i] = in[i]*exp(ii*m_omega[i]*m_mint)/m_L;
    // Stuck with kiss_fft_cpx definition of complex -> use reinterpret and assume 
    // kiss interface won't change :(
    kiss_fft(m_cfg_reverse,reinterpret_cast<const kiss_fft_cpx*>(m_temp.data()),\
                  reinterpret_cast<kiss_fft_cpx*>(out));
}

}