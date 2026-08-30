#pragma once

#include "spida/grid/uniformCVT.h"
#include "spida/helper/constants.h"

#include <vector>

// kiss_fft_cfg is kissfft's own opaque pointer typedef (struct
// kiss_fft_state*, never defined even by kiss_fft.h itself) — forward
// declared here so this public header doesn't leak a kissfft header into
// every consumer's translation unit. fftCVT.cpp includes the real
// kiss_fft.h for the definition it actually needs.
struct kiss_fft_state;
using kiss_fft_cfg = kiss_fft_state*;

namespace spida {

class FFTCVT {
public:
    explicit FFTCVT(const UniformGridCVT& grid);
    ~FFTCVT();
    FFTCVT() = delete;
    FFTCVT(const FFTCVT& sp) = delete;
    FFTCVT& operator=(const FFTCVT& sp) = delete;
    void T_To_ST(const dcmplx* in, dcmplx* out) noexcept;
    void ST_To_T(const dcmplx* in, dcmplx* out) noexcept;

    void T_To_ST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) noexcept
    {
        T_To_ST(in.data(), out.data());
    }

    void ST_To_T(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) noexcept
    {
        ST_To_T(in.data(), out.data());
    }

private:
    unsigned m_nt;
    double m_mint;
    double m_L;
    std::vector<double> m_omega;
    std::vector<dcmplx> m_temp;

    kiss_fft_cfg m_cfg_forward;
    kiss_fft_cfg m_cfg_reverse;
};

} // namespace spida