#pragma once

#include "spida/grid/uniformCVT.h"
#include "spida/helper/constants.h"

#include <cassert>
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
    // Move is suppressed too, not just implicitly left unavailable by the
    // deleted copy ops above -- m_cfg_forward/m_cfg_reverse are raw owning
    // handles freed in the destructor, and a moved-from instance's
    // destructor would double-free them without an explicit transfer this
    // class doesn't implement.
    FFTCVT(FFTCVT&& sp) = delete;
    FFTCVT& operator=(FFTCVT&& sp) = delete;
    void T_To_ST(const dcmplx* in, dcmplx* out) noexcept;
    void ST_To_T(const dcmplx* in, dcmplx* out) noexcept;

    void T_To_ST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) noexcept
    {
        // Debug-only guard against caller/grid size drift -- the
        // raw-pointer overload this delegates to has no length to check
        // itself (no size parameter at all; m_nt is trusted implicitly).
        // Compiled out under NDEBUG, so no Release-path cost. ST domain is
        // the same size as T for a complex-to-complex transform (CVT's
        // getNst() == getNt()).
        assert(in.size() == m_nt);
        assert(out.size() == m_nt);
        T_To_ST(in.data(), out.data());
    }

    void ST_To_T(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) noexcept
    {
        assert(in.size() == m_nt);
        assert(out.size() == m_nt);
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