#pragma once

// Real-valued KdV (u_t + 6 u u_x + u_xxx = 0), promoted from the spida-worker
// implementation originally authored in spida-console (see
// docs/adr/0003-worker-relocation-and-cooperative-cancellation.md) — verified
// there against the exact single-soliton closed-form solution (max abs error
// 0.0024 at n=1024, consistent with grid spacing, not drift). Mirrors
// spida/models/burgers.h's structure exactly: a model class owning L()/NL(),
// paired with a PropagatorCV reporting "X"/"SX".
//
// Also holds KdvCv/KdvCvPropagator below (ModelKind::kdv_cv) — the exact
// same PDE, promoted from demos/kdvCV.cpp, solved on a full-complex-FFT
// uniform grid (GridKind::uniform_cvx, SpidaCVX) instead of Kdv's
// real-optimized half-spectrum one (SpidaRVX). Kept in this file rather
// than a separate one since it's the same equation family, not a
// different model — see docs/adr/0001-spida-console-backend-groundwork.md's
// Phase D addendum for why kdv_cv exists as its own ModelKind at all
// (the proposal's own domain.ts lists it separately from kdv_rv).

#include "spida/CVX.h"
#include "spida/RVX.h"
#include "spida/grid/uniformCVX.h"
#include "spida/grid/uniformRVX.h"
#include "spida/helper/constants.h"
#include "spida/propagator/propagator.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

#include <utils/report.hpp>

namespace spida::models {

/// u_t + 6 u u_x + u_xxx = 0, periodic on [a, b] — the standard
/// normalization, no free PDE coefficient (unlike Burgers' mu). Solved the
/// same pseudospectral way:
///   L(k)  = i k^3    (dispersive: -u_xxx in Fourier space, (ik)^3 = -i k^3)
///   NL(u) = -6 u u_x
/// L is purely imaginary (dispersive, not dissipative like Burgers'
/// -mu*k^2) — ETD35 makes no assumption either way, its LinOp is a
/// std::vector<dcmplx> throughout.
class Kdv {
public:
    explicit Kdv(const spida::UniformGridRVX& grid)
        : m_grid(grid),
          m_spi(grid),
          m_uphys(grid.getNx()),
          m_uxphys(grid.getNx()),
          m_uxsp(grid.getNsx()),
          m_L(grid.getNsx())
    {
        const auto& sx = grid.getSX();
        for (std::size_t i = 0; i < sx.size(); i++)
            m_L[i] = spida::dcmplx(0.0, std::pow(sx[i], 3));
        m_NL = [this](const std::vector<spida::dcmplx>& in, std::vector<spida::dcmplx>& out) {
            m_spi.SX_To_X(in, m_uphys);
            m_spi.dSX(in, m_uxsp);
            m_spi.SX_To_X(m_uxsp, m_uxphys);
            for (std::size_t i = 0; i < m_grid.getNx(); i++)
                m_uphys[i] = -6.0 * m_uphys[i] * m_uxphys[i];
            m_spi.X_To_SX(m_uphys, out);
        };
    }

    [[nodiscard]] std::vector<spida::dcmplx>& L()
    {
        return m_L;
    }

    [[nodiscard]] std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)>&
    NL()
    {
        return m_NL;
    }

    [[nodiscard]] spida::SpidaRVX& spida()
    {
        return m_spi;
    }

private:
    spida::UniformGridRVX m_grid;
    spida::SpidaRVX m_spi;
    std::vector<double> m_uphys;
    std::vector<double> m_uxphys;
    std::vector<spida::dcmplx> m_uxsp;
    std::vector<spida::dcmplx> m_L;
    std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)> m_NL;
};

/// Single-soliton initial condition: u(x,0) = (c/2) sech^2(sqrt(c)/2 (x-x0)),
/// x0 = grid midpoint — the exact solution translates at speed c with no
/// change of shape, u(x,t) = (c/2) sech^2(sqrt(c)/2 (x - c*t - x0)). `c`
/// (solitonSpeed) is the one free parameter, since KdV's PDE itself has none.
class KdvPropagator : public spida::PropagatorCV {
public:
    KdvPropagator(const std::filesystem::path& path, Kdv& model, double solitonSpeed)
        : PropagatorCV(path),
          m_spi(model.spida()),
          m_usp(model.spida().getGridX().getNsx(), 0.0),
          m_uphys(model.spida().getGridX().getNx(), 0.0)
    {
        const auto& x = m_spi.getX();
        const double x0 = (x.front() + x.back()) / 2.0;
        const double amp = solitonSpeed / 2.0;
        const double width = std::sqrt(solitonSpeed) / 2.0;
        for (std::size_t i = 0; i < x.size(); i++) {
            const double s = 1.0 / std::cosh(width * (x[i] - x0));
            m_uphys[i] = amp * s * s;
        }
        m_spi.X_To_SX(m_uphys, m_usp);
        initReport();
    }

    ~KdvPropagator() override = default;

    void updateFields(double) override
    {
        m_spi.SX_To_X(m_usp, m_uphys);
    }

    [[nodiscard]] std::vector<spida::dcmplx>& propagator() override
    {
        return m_usp;
    }

    bool checkDiverged(double) override
    {
        for (double v : m_uphys) {
            if (!std::isfinite(v))
                return true;
        }
        return false;
    }

private:
    void initReport()
    {
        const auto& x = m_spi.getGridX().getX();
        auto report = std::make_unique<spida::Report1D<double, double>>("X", x, m_uphys);
        PropagatorCV::addReport(std::move(report));
        const auto& sx = m_spi.getGridX().getSX();
        auto reportsp =
            std::make_unique<spida::ReportComplex1D<double, double>>("SX", sx, m_usp);
        PropagatorCV::addReport(std::move(reportsp));
    }

    spida::SpidaRVX& m_spi;
    std::vector<spida::dcmplx> m_usp;
    std::vector<double> m_uphys;
};

/// Same PDE as Kdv above, solved via SpidaCVX's full complex FFT on a
/// GridKind::uniform_cvx grid instead of SpidaRVX's real-optimized
/// half-spectrum one. L(k) = i k^3, same formula as Kdv's — the only
/// difference is which transform computes it over.
class KdvCv {
public:
    explicit KdvCv(const spida::UniformGridCVX& grid)
        : m_grid(grid),
          m_spi(grid),
          m_uphys(grid.getNx()),
          m_uxphys(grid.getNx()),
          m_uxsp(grid.getNsx()),
          m_L(grid.getNsx())
    {
        const auto& sx = grid.getSX();
        for (std::size_t i = 0; i < sx.size(); i++)
            m_L[i] = spida::dcmplx(0.0, std::pow(sx[i], 3));
        m_NL = [this](const std::vector<spida::dcmplx>& in, std::vector<spida::dcmplx>& out) {
            m_spi.SX_To_X(in, m_uphys);
            m_spi.dSX(in, m_uxsp);
            m_spi.SX_To_X(m_uxsp, m_uxphys);
            for (std::size_t i = 0; i < m_grid.getNx(); i++)
                m_uphys[i] = -6.0 * m_uphys[i] * m_uxphys[i];
            m_spi.X_To_SX(m_uphys, out);
        };
    }

    [[nodiscard]] std::vector<spida::dcmplx>& L()
    {
        return m_L;
    }

    [[nodiscard]] std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)>&
    NL()
    {
        return m_NL;
    }

    [[nodiscard]] spida::SpidaCVX& spida()
    {
        return m_spi;
    }

private:
    spida::UniformGridCVX m_grid;
    spida::SpidaCVX m_spi;
    std::vector<spida::dcmplx> m_uphys;
    std::vector<spida::dcmplx> m_uxphys;
    std::vector<spida::dcmplx> m_uxsp;
    std::vector<spida::dcmplx> m_L;
    std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)> m_NL;
};

/// Fixed 5-soliton initial condition, matching demos/kdvCV.cpp exactly —
/// no modelParams (same precedent as Ks having none): this model exists to
/// demonstrate the complex-valued-grid pipeline on a nontrivial
/// multi-soliton field, not as a user-tunable scenario. Soliton centers
/// span x in [-120, 0], so this model expects a domain at least that wide
/// (the demo uses [-150, 150]) — same "needs room" caveat KdvPropagator's
/// own header comment gives for its single soliton.
///
/// Numerical verification: since both the PDE (real coefficients) and this
/// initial condition are real-valued in physical space, u(x,t) should stay
/// real for all t under exact evolution — a standard fact about real-
/// coefficient PDEs. propagator() exposes the SPECTRAL array though, which
/// is generally complex at nonzero frequencies even for a real signal —
/// the actual checkable invariant is that array's Hermitian symmetry,
/// usp[N-k] == conj(usp[k]). Checked for real against the actual solver
/// output at two points (right after construction, and after run()) — see
/// test/config_tests.cpp's KDV_CV_STAYS_REAL_VALUED for the measured
/// values (~5e-16 relative at t0, confirming the IC/transform introduce no
/// asymmetry; ~2.6e-6 at tf=0.01, a legitimate time-stepping accumulation).
class KdvCvPropagator : public spida::PropagatorCV {
public:
    KdvCvPropagator(const std::filesystem::path& path, KdvCv& model)
        : PropagatorCV(path),
          m_spi(model.spida()),
          m_usp(model.spida().getGridX().getNsx(), 0.0),
          m_uphys(model.spida().getGridX().getNx(), 0.0)
    {
        static const std::vector<double> A0{0.6, 0.5, 0.4, 0.3, 0.2};
        static const std::vector<double> x0{-120.0, -90.0, -60.0, -30.0, 0.0};
        const auto& x = m_spi.getX();
        for (std::size_t i = 0; i < x.size(); i++) {
            for (std::size_t k = 0; k < A0.size(); k++) {
                const double s = 1.0 / std::cosh(A0[k] * (x[i] - x0[k]) / 2.0);
                m_uphys[i] += 0.5 * A0[k] * A0[k] * s * s;
            }
        }
        m_spi.X_To_SX(m_uphys, m_usp);
        initReport();
    }

    ~KdvCvPropagator() override = default;

    void updateFields(double) override
    {
        m_spi.SX_To_X(m_usp, m_uphys);
    }

    [[nodiscard]] std::vector<spida::dcmplx>& propagator() override
    {
        return m_usp;
    }

    bool checkDiverged(double) override
    {
        for (const auto& v : m_uphys) {
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
                return true;
        }
        return false;
    }

private:
    void initReport()
    {
        const auto& x = m_spi.getGridX().getX();
        auto report = std::make_unique<spida::ReportComplex1D<double, double>>("X", x, m_uphys);
        PropagatorCV::addReport(std::move(report));
        const auto& sx = m_spi.getGridX().getSX();
        auto reportsp =
            std::make_unique<spida::ReportComplex1D<double, double>>("SX", sx, m_usp);
        PropagatorCV::addReport(std::move(reportsp));
    }

    spida::SpidaCVX& m_spi;
    std::vector<spida::dcmplx> m_usp;
    std::vector<spida::dcmplx> m_uphys;
};

} // namespace spida::models
