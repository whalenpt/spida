
#include <cmath>
#include <random>
#include <vector>

#include <gtest/gtest.h>
#include <spida/CVT.h>
#include <spida/CVX.h>
#include <spida/grid/uniformCVT.h>
#include <spida/grid/uniformCVX.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>
#include <spida/RVT.h>
#include <spida/RVX.h>
#include <utils/except.h>
#include <utils/math.hpp>

using spida::dcmplx;
using spida::ii;
using spida::PI;

// ============================================================
// SpidaCVX wrapper tests
// ============================================================

TEST(SPIDA_CVX_TEST, ACCESSORS)
{
    unsigned N = 32;
    spida::UniformGridCVX grid(N, -3.0, 3.0);
    spida::SpidaCVX spi(grid);

    EXPECT_EQ(spi.getGridX().getNx(), N);
    EXPECT_DOUBLE_EQ(spi.getGridX().getMinX(), -3.0);
    EXPECT_DOUBLE_EQ(spi.getGridX().getMaxX(), 3.0);
    EXPECT_EQ(spi.getX().size(), N);
    EXPECT_EQ(spi.getSX().size(), N);
    EXPECT_EQ(spi.getX(), grid.getX());
    EXPECT_EQ(spi.getSX(), grid.getSX());
}

TEST(SPIDA_CVX_TEST, TRANSFORM_ACCESSOR)
{
    unsigned N = 16;
    spida::UniformGridCVX grid(N, -1.0, 1.0);
    spida::SpidaCVX spi(grid);
    std::vector<dcmplx> in(N, 1.0), sp(N), out(N);
    spi.getTransformX().X_To_SX(in, sp);
    spi.getTransformX().SX_To_X(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_CVX_TEST, X_TO_SX_AND_BACK)
{
    unsigned N = 64;
    spida::UniformGridCVX grid(N, -6.0, 6.0);
    spida::SpidaCVX spi(grid);

    std::default_random_engine gen(42);
    std::normal_distribution<double> dist;
    std::vector<dcmplx> in(N), sp(N), out(N);
    for (auto& v : in)
        v = dcmplx(dist(gen), dist(gen));

    spi.X_To_SX(in, sp);
    spi.SX_To_X(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_CVX_TEST, X_TO_SX_GAUSS)
{
    // F{exp(-a*x^2)} = sqrt(pi/a) * exp(-kx^2 / (4a))
    unsigned N = 64;
    double a = 2.0;
    spida::UniformGridCVX grid(N, -6.0, 6.0);
    spida::SpidaCVX spi(grid);

    const auto& x = spi.getX();
    const auto& kx = spi.getSX();
    std::vector<dcmplx> in(N), sp(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::exp(-a * x[i] * x[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = std::sqrt(PI / a) * std::exp(-kx[i] * kx[i] / (4.0 * a));

    spi.X_To_SX(in, sp);
    EXPECT_LT(detail::relative_error(expect, sp), 1e-5);
}

TEST(SPIDA_CVX_TEST, DX_NOOP_N0)
{
    unsigned N = 16;
    spida::UniformGridCVX grid(N, -1.0, 1.0);
    spida::SpidaCVX spi(grid);

    std::vector<dcmplx> in(N, dcmplx(1.0, 0.5));
    std::vector<dcmplx> out(N, dcmplx(99.0, 99.0));
    spi.dX(in, out, 0);
    // n=0 should be a no-op; out must be unchanged
    for (unsigned i = 0; i < N; i++)
        EXPECT_EQ(out[i], dcmplx(99.0, 99.0));
}

TEST(SPIDA_CVX_TEST, DSX_NOOP_N0)
{
    unsigned N = 16;
    spida::UniformGridCVX grid(N, -1.0, 1.0);
    spida::SpidaCVX spi(grid);

    std::vector<dcmplx> in(N, dcmplx(1.0, 0.5));
    std::vector<dcmplx> out(N, dcmplx(99.0, 99.0));
    spi.dSX(in, out, 0);
    for (unsigned i = 0; i < N; i++)
        EXPECT_EQ(out[i], dcmplx(99.0, 99.0));
}

TEST(SPIDA_CVX_TEST, DX_FIRST_DERIVATIVE_SIN)
{
    // d/dx sin(x) = cos(x)
    unsigned N = 64;
    spida::UniformGridCVX grid(N, 0.0, 4.0 * PI);
    spida::SpidaCVX spi(grid);

    const auto& x = spi.getX();
    std::vector<dcmplx> in(N), out(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(x[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = std::cos(x[i]);

    spi.dX(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}

TEST(SPIDA_CVX_TEST, DX_SECOND_DERIVATIVE_SIN)
{
    // d^2/dx^2 sin(x) = -sin(x)
    unsigned N = 64;
    spida::UniformGridCVX grid(N, 0.0, 4.0 * PI);
    spida::SpidaCVX spi(grid);

    const auto& x = spi.getX();
    std::vector<dcmplx> in(N), out(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(x[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = -std::sin(x[i]);

    spi.dX(in, out, 2);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}

TEST(SPIDA_CVX_TEST, DSX_MATCHES_DX_SPECTRAL_STEP)
{
    // dSX applied to F{sin(x)} should equal F{cos(x)}
    unsigned N = 64;
    spida::UniformGridCVX grid(N, 0.0, 4.0 * PI);
    spida::SpidaCVX spi(grid);

    const auto& x = spi.getX();
    std::vector<dcmplx> in_sp(N), out_sp(N), expect_sp(N);

    // Compute F{sin(x)} and F{cos(x)} to compare
    std::vector<dcmplx> sin_x(N), cos_x(N);
    for (unsigned i = 0; i < N; i++)
        sin_x[i] = std::sin(x[i]);
    for (unsigned i = 0; i < N; i++)
        cos_x[i] = std::cos(x[i]);
    spi.X_To_SX(sin_x, in_sp);
    spi.X_To_SX(cos_x, expect_sp);

    spi.dSX(in_sp, out_sp);
    EXPECT_LT(detail::relative_error(expect_sp, out_sp), 1e-5);
}

// ============================================================
// SpidaRVX wrapper tests
// ============================================================

TEST(SPIDA_RVX_TEST, ACCESSORS)
{
    unsigned N = 32;
    spida::UniformGridRVX grid(N, -3.0, 3.0);
    spida::SpidaRVX spi(grid);

    EXPECT_EQ(spi.getGridX().getNx(), N);
    EXPECT_DOUBLE_EQ(spi.getGridX().getMinX(), -3.0);
    EXPECT_DOUBLE_EQ(spi.getGridX().getMaxX(), 3.0);
    EXPECT_EQ(spi.getX().size(), N);
    EXPECT_EQ(spi.getSX().size(), grid.getNsx());
    EXPECT_EQ(spi.getX(), grid.getX());
    EXPECT_EQ(spi.getSX(), grid.getSX());
}

TEST(SPIDA_RVX_TEST, TRANSFORM_ACCESSOR)
{
    unsigned N = 16;
    spida::UniformGridRVX grid(N, -1.0, 1.0);
    spida::SpidaRVX spi(grid);
    std::vector<double> in(N, 1.0), out(N);
    std::vector<dcmplx> sp(grid.getNsx());
    spi.getTransformX().X_To_SX(in, sp);
    spi.getTransformX().SX_To_X(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_RVX_TEST, X_TO_SX_AND_BACK)
{
    unsigned N = 64;
    spida::UniformGridRVX grid(N, -6.0, 6.0);
    spida::SpidaRVX spi(grid);

    std::default_random_engine gen(7);
    std::normal_distribution<double> dist;
    std::vector<double> in(N), out(N);
    std::vector<dcmplx> sp(grid.getNsx());
    for (auto& v : in)
        v = dist(gen);

    spi.X_To_SX(in, sp);
    spi.SX_To_X(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_RVX_TEST, X_TO_SX_GAUSS)
{
    // F{exp(-a*x^2)} = sqrt(pi/a) * exp(-kx^2 / (4a))
    unsigned N = 64;
    double a = 2.0;
    spida::UniformGridRVX grid(N, -6.0, 6.0);
    spida::SpidaRVX spi(grid);

    const auto& x = spi.getX();
    const auto& kx = spi.getSX();
    std::vector<double> in(N);
    std::vector<dcmplx> sp(grid.getNsx()), expect(grid.getNsx());
    for (unsigned i = 0; i < N; i++)
        in[i] = std::exp(-a * x[i] * x[i]);
    for (unsigned i = 0; i < kx.size(); i++)
        expect[i] = std::sqrt(PI / a) * std::exp(-kx[i] * kx[i] / (4.0 * a));

    spi.X_To_SX(in, sp);
    EXPECT_LT(detail::relative_error(expect, sp), 1e-5);
}

TEST(SPIDA_RVX_TEST, DX_NOOP_N0)
{
    unsigned N = 16;
    spida::UniformGridRVX grid(N, -1.0, 1.0);
    spida::SpidaRVX spi(grid);

    std::vector<double> in(N, 3.14), out(N, 99.0);
    spi.dX(in, out, 0);
    for (unsigned i = 0; i < N; i++)
        EXPECT_DOUBLE_EQ(out[i], 99.0);
}

TEST(SPIDA_RVX_TEST, DSX_TWO_BUFFER_OVERLOAD)
{
    // dSX(in, out, n) — two-buffer version
    unsigned N = 64;
    spida::UniformGridRVX grid(N, 0.0, 4.0 * PI);
    spida::SpidaRVX spi(grid);

    const auto& x = spi.getX();
    std::vector<double> sin_x(N), cos_x(N);
    for (unsigned i = 0; i < N; i++)
        sin_x[i] = std::sin(x[i]);
    for (unsigned i = 0; i < N; i++)
        cos_x[i] = std::cos(x[i]);

    std::vector<dcmplx> in_sp(grid.getNsx()), out_sp(grid.getNsx()), expect_sp(grid.getNsx());
    spi.X_To_SX(sin_x, in_sp);
    spi.X_To_SX(cos_x, expect_sp);

    spi.dSX(in_sp, out_sp);
    EXPECT_LT(detail::relative_error(expect_sp, out_sp), 1e-5);
}

TEST(SPIDA_RVX_TEST, DSX_INPLACE_NOOP_N0)
{
    unsigned N = 16;
    spida::UniformGridRVX grid(N, -1.0, 1.0);
    spida::SpidaRVX spi(grid);

    std::vector<dcmplx> v(grid.getNsx(), dcmplx(1.5, 2.5));
    spi.dSX(v, 0);
    for (const auto& c : v)
        EXPECT_EQ(c, dcmplx(1.5, 2.5));
}

TEST(SPIDA_RVX_TEST, DSX_INPLACE_MATCHES_TWO_BUFFER)
{
    unsigned N = 32;
    spida::UniformGridRVX grid(N, 0.0, 2.0 * PI);
    spida::SpidaRVX spi(grid);

    const auto& x = spi.getX();
    std::vector<double> in(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(x[i]);

    std::vector<dcmplx> sp(grid.getNsx());
    spi.X_To_SX(in, sp);

    std::vector<dcmplx> out_two(grid.getNsx()), out_inplace = sp;
    spi.dSX(sp, out_two);
    spi.dSX(out_inplace);

    EXPECT_LT(detail::relative_error(out_two, out_inplace), 1e-12);
}

// ============================================================
// SpidaCVT wrapper tests
// ============================================================

TEST(SPIDA_CVT_TEST, ACCESSORS)
{
    unsigned N = 32;
    spida::UniformGridCVT grid(N, -2.0, 2.0);
    spida::SpidaCVT spi(grid);

    EXPECT_EQ(spi.getGridT().getNt(), N);
    EXPECT_DOUBLE_EQ(spi.getGridT().getMinT(), -2.0);
    EXPECT_DOUBLE_EQ(spi.getGridT().getMaxT(), 2.0);
    EXPECT_EQ(spi.getT().size(), N);
    EXPECT_EQ(spi.getST().size(), N);
    EXPECT_EQ(spi.getT(), grid.getT());
    EXPECT_EQ(spi.getST(), grid.getST());
}

TEST(SPIDA_CVT_TEST, TRANSFORM_ACCESSOR)
{
    unsigned N = 16;
    spida::UniformGridCVT grid(N, -1.0, 1.0);
    spida::SpidaCVT spi(grid);
    std::vector<dcmplx> in(N, 1.0), sp(N), out(N);
    spi.getTransformT().T_To_ST(in, sp);
    spi.getTransformT().ST_To_T(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_CVT_TEST, T_TO_ST_AND_BACK)
{
    unsigned N = 64;
    spida::UniformGridCVT grid(N, -6.0, 6.0);
    spida::SpidaCVT spi(grid);

    std::default_random_engine gen(13);
    std::normal_distribution<double> dist;
    std::vector<dcmplx> in(N), sp(N), out(N);
    for (auto& v : in)
        v = dcmplx(dist(gen), dist(gen));

    spi.T_To_ST(in, sp);
    spi.ST_To_T(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_CVT_TEST, T_TO_ST_GAUSS)
{
    // F{exp(-a*t^2)} = sqrt(pi/a) * exp(-omega^2 / (4a))  [with exp(+i*omega*t) convention]
    unsigned N = 64;
    double a = 2.0;
    spida::UniformGridCVT grid(N, -6.0, 6.0);
    spida::SpidaCVT spi(grid);

    const auto& t = spi.getT();
    const auto& omega = spi.getST();
    std::vector<dcmplx> in(N), sp(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::exp(-a * t[i] * t[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = std::sqrt(PI / a) * std::exp(-omega[i] * omega[i] / (4.0 * a));

    spi.T_To_ST(in, sp);
    EXPECT_LT(detail::relative_error(expect, sp), 1e-5);
}

TEST(SPIDA_CVT_TEST, DT_NOOP_N0)
{
    unsigned N = 16;
    spida::UniformGridCVT grid(N, -1.0, 1.0);
    spida::SpidaCVT spi(grid);

    std::vector<dcmplx> in(N, dcmplx(2.0, -1.0));
    std::vector<dcmplx> out(N, dcmplx(99.0, 99.0));
    spi.dT(in, out, 0);
    for (unsigned i = 0; i < N; i++)
        EXPECT_EQ(out[i], dcmplx(99.0, 99.0));
}

TEST(SPIDA_CVT_TEST, DST_NOOP_N0)
{
    unsigned N = 16;
    spida::UniformGridCVT grid(N, -1.0, 1.0);
    spida::SpidaCVT spi(grid);

    std::vector<dcmplx> in(N, dcmplx(2.0, -1.0));
    std::vector<dcmplx> out(N, dcmplx(99.0, 99.0));
    spi.dST(in, out, 0);
    for (unsigned i = 0; i < N; i++)
        EXPECT_EQ(out[i], dcmplx(99.0, 99.0));
}

TEST(SPIDA_CVT_TEST, DT_FIRST_DERIVATIVE_SIN)
{
    // d/dt sin(t) = cos(t)
    unsigned N = 64;
    spida::UniformGridCVT grid(N, 0.0, 4.0 * PI);
    spida::SpidaCVT spi(grid);

    const auto& t = spi.getT();
    std::vector<dcmplx> in(N), out(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(t[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = std::cos(t[i]);

    spi.dT(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}

TEST(SPIDA_CVT_TEST, DT_SECOND_DERIVATIVE_SIN)
{
    // d^2/dt^2 sin(t) = -sin(t)
    unsigned N = 64;
    spida::UniformGridCVT grid(N, 0.0, 4.0 * PI);
    spida::SpidaCVT spi(grid);

    const auto& t = spi.getT();
    std::vector<dcmplx> in(N), out(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(t[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = -std::sin(t[i]);

    spi.dT(in, out, 2);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}

TEST(SPIDA_CVT_TEST, DST_MATCHES_DT_SPECTRAL_STEP)
{
    // dST of F{sin(t)} should equal F{cos(t)}
    unsigned N = 64;
    spida::UniformGridCVT grid(N, 0.0, 4.0 * PI);
    spida::SpidaCVT spi(grid);

    const auto& t = spi.getT();
    std::vector<dcmplx> sin_t(N), cos_t(N);
    for (unsigned i = 0; i < N; i++)
        sin_t[i] = std::sin(t[i]);
    for (unsigned i = 0; i < N; i++)
        cos_t[i] = std::cos(t[i]);

    std::vector<dcmplx> in_sp(N), out_sp(N), expect_sp(N);
    spi.T_To_ST(sin_t, in_sp);
    spi.T_To_ST(cos_t, expect_sp);

    spi.dST(in_sp, out_sp);
    EXPECT_LT(detail::relative_error(expect_sp, out_sp), 1e-5);
}

// ============================================================
// SpidaRVT wrapper tests
// ============================================================

TEST(SPIDA_RVT_TEST, ACCESSORS)
{
    unsigned N = 64;
    double minT = -3.0, maxT = 3.0;
    spida::UniformGridRVT grid(N, minT, maxT);
    spida::SpidaRVT spi(grid);

    EXPECT_EQ(spi.getGridT().getNt(), N);
    EXPECT_DOUBLE_EQ(spi.getGridT().getMinT(), minT);
    EXPECT_DOUBLE_EQ(spi.getGridT().getMaxT(), maxT);
    EXPECT_EQ(spi.getT().size(), N);
    EXPECT_EQ(spi.getST().size(), grid.getNst());
    EXPECT_EQ(spi.getT(), grid.getT());
    EXPECT_EQ(spi.getST(), grid.getST());
}

TEST(SPIDA_RVT_TEST, TRANSFORM_ACCESSOR)
{
    unsigned N = 32;
    spida::UniformGridRVT grid(N, -1.0, 1.0);
    spida::SpidaRVT spi(grid);
    std::vector<double> in(N, 1.0), out(N);
    std::vector<dcmplx> sp(grid.getNst());
    spi.getTransformT().T_To_ST(in, sp);
    spi.getTransformT().ST_To_T(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_RVT_TEST, T_TO_ST_AND_BACK)
{
    unsigned N = 64;
    spida::UniformGridRVT grid(N, -6.0, 6.0);
    spida::SpidaRVT spi(grid);

    std::default_random_engine gen(3);
    std::normal_distribution<double> dist;
    std::vector<double> in(N), out(N);
    std::vector<dcmplx> sp(grid.getNst());
    for (auto& v : in)
        v = dist(gen);

    spi.T_To_ST(in, sp);
    spi.ST_To_T(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(SPIDA_RVT_TEST, CVT_TO_ST_AND_BACK)
{
    // CVT_To_ST keeps only positive frequencies, so the round-trip is only
    // lossless for band-limited (analytic) signals.  Build such a signal by
    // going spectral -> time first, then verify time -> spectral recovers it.
    unsigned N = 64;
    spida::UniformGridRVT grid(N, -6.0, 6.0);
    spida::SpidaRVT spi(grid);

    unsigned nst = grid.getNst();
    std::vector<dcmplx> sp_in(nst, 0.0);
    sp_in[nst / 2] = 1.0; // single mid-band frequency: guaranteed analytic

    std::vector<dcmplx> time_domain(N), sp_out(nst);
    spi.ST_To_CVT(sp_in, time_domain);
    spi.CVT_To_ST(time_domain, sp_out);
    EXPECT_LT(detail::relative_error(sp_in, sp_out), 1e-10);
}

TEST(SPIDA_RVT_TEST, T_TO_ST_GAUSS)
{
    // F{exp(-a*t^2)} = sqrt(pi/a) * exp(-omega^2 / (4a))
    unsigned N = 64;
    double a = 2.0;
    spida::UniformGridRVT grid(N, -6.0, 6.0);
    spida::SpidaRVT spi(grid);

    const auto& t = spi.getT();
    const auto& omega = spi.getST();
    std::vector<double> in(N);
    std::vector<dcmplx> sp(grid.getNst()), expect(grid.getNst());
    for (unsigned i = 0; i < N; i++)
        in[i] = std::exp(-a * t[i] * t[i]);
    for (unsigned i = 0; i < omega.size(); i++)
        expect[i] = std::sqrt(PI / a) * std::exp(-omega[i] * omega[i] / (4.0 * a));

    spi.T_To_ST(in, sp);
    EXPECT_LT(detail::relative_error(expect, sp), 1e-5);
}

TEST(SPIDA_RVT_TEST, DT_NEGATIVE_N_THROWS)
{
    unsigned N = 32;
    spida::UniformGridRVT grid(N, 0.0, 2.0 * PI);
    spida::SpidaRVT spi(grid);

    std::vector<double> in(N, 1.0), out(N);
    EXPECT_THROW(spi.dT(in, out, -1), detail::Exception);
}

TEST(SPIDA_RVT_TEST, DT_FIRST_DERIVATIVE_SIN)
{
    // d/dt sin(t) = cos(t)
    unsigned N = 64;
    double tmin = 0.0, tmax = 2.0 * PI;
    spida::UniformGridRVT grid(N, tmin, tmax);
    spida::SpidaRVT spi(grid);

    const auto& t = spi.getT();
    std::vector<double> in(N), out(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(t[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = std::cos(t[i]);

    spi.dT(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}

TEST(SPIDA_RVT_TEST, DT_SECOND_DERIVATIVE_SIN)
{
    // d^2/dt^2 sin(t) = -sin(t)
    unsigned N = 64;
    spida::UniformGridRVT grid(N, 0.0, 2.0 * PI);
    spida::SpidaRVT spi(grid);

    const auto& t = spi.getT();
    std::vector<double> in(N), out(N), expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = std::sin(t[i]);
    for (unsigned i = 0; i < N; i++)
        expect[i] = -std::sin(t[i]);

    spi.dT(in, out, 2);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}
