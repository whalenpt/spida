#pragma once

// Real-valued KdV (u_t + 6 u u_x + u_xxx = 0), promoted from the spida-worker
// implementation originally authored in spida-console (see
// docs/adr/0003-worker-relocation-and-cooperative-cancellation.md) — verified
// there against the exact single-soliton closed-form solution (max abs error
// 0.0024 at n=1024, consistent with grid spacing, not drift). Mirrors
// spida/models/burgers.h's structure exactly: a model class owning L()/NL(),
// paired with a PropagatorCV reporting "X"/"SX".

#include "spida/RVX.h"
#include "spida/grid/uniformRVX.h"
#include "spida/helper/constants.h"
#include "spida/propagator/propagator.h"

#include <cmath>
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

} // namespace spida::models
