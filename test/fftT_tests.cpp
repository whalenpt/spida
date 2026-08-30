

#include <fstream>
#include <random>
#include <stdexcept>

#include <gtest/gtest.h>
#include <spida/CVT.h>
#include <spida/grid/uniformCVT.h>
#include <spida/grid/uniformRVT.h>
#include <spida/helper/constants.h>
#include <spida/RVT.h>
#include <spida/shape/shapeT.h>
#include <spida/transform/fftCVT.h>
#include <spida/transform/fftRVT.h>
#include <utils/math.hpp>

// FFTCVT defined such that F{f(t)} = \integral_{-\inf}^{\inf}f(t)exp(i*omega*t) dt
// Test that forward fft followed by inverse fft yields identity
TEST(FFTCVT_TEST, INVERSES)
{
    unsigned N = 32;
    using spida::dcmplx;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0, 1.0);
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::FFTCVT tr(spida::UniformGridCVT{N, -1, 1});
    tr.T_To_ST(in, out);
    tr.ST_To_T(out, expect);

    EXPECT_LT(detail::relative_error(in, expect), 1e-10);
}

// a = 1/tp^2
// F{exp(-a*t^2)} = sqrt(pi/a)*exp(-omega^2/(4*a))
// F{exp(-(t/tp)^2)}= tp*sqrt(pi)*exp(-tp^2*omega^2/4)
TEST(FFTCVT_TEST, GAUSS)
{
    unsigned N = 64;
    using spida::dcmplx;
    using spida::ii;
    using spida::PI;

    std::vector<dcmplx> in(N, 0.0);
    std::vector<dcmplx> out(N, 0.0);
    std::vector<dcmplx> expect(N, 0.0);

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridCVT grid(N, xmin, xmax);
    const std::vector<double> t = grid.getT();
    const std::vector<double> omega = grid.getST();

    double a = 2.0;
    for (size_t i = 0; i < t.size(); i++)
        in[i] = exp(-a * pow(t[i], 2));
    for (size_t i = 0; i < omega.size(); i++)
        expect[i] = sqrt(PI / a) * exp(-pow(omega[i], 2) / (4.0 * a));

    spida::FFTCVT tr(grid);
    tr.T_To_ST(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-5);
}

TEST(FFTCVT_TEST, COS)
{
    // F{cos(at)} = PI*(delta(omega-a) + delta(omega+a))
    unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    spida::UniformGridCVT grid(N, 0.0, 2.0 * PI);
    const std::vector<double> t = grid.getT();
    for (size_t i = 0; i < t.size(); i++)
        in[i] = cos(8 * t[i]);

    spida::FFTCVT tr(grid);
    tr.T_To_ST(in, out);
    EXPECT_DOUBLE_EQ(out[8].real(), PI);
}

TEST(FFTCVT_TEST, DERIVATIVE_SIN)
{
    unsigned N = 32;

    using spida::dcmplx;
    using spida::PI;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVT grid(N, 0, 2 * PI);
    const std::vector<double> t = grid.getT();
    for (size_t i = 0; i < t.size(); i++)
        in[i] = sin(t[i]);
    for (size_t i = 0; i < t.size(); i++)
        expect[i] = cos(t[i]);

    spida::SpidaCVT spi{grid};
    spi.dT(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-6);
}

TEST(FFTCVT_TEST, DERIVATIVE_GAUSS)
{
    unsigned N = 32;

    using spida::dcmplx;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVT grid(N, -6, 6);
    const std::vector<double> t = grid.getT();
    for (size_t i = 0; i < t.size(); i++)
        in[i] = exp(-pow(t[i], 2));
    for (size_t i = 0; i < t.size(); i++)
        expect[i] = -2.0 * t[i] * exp(-pow(t[i], 2));

    spida::SpidaCVT spi{grid};
    spi.dT(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-6);
}

// FFTRVT defined such that F{f(t)} = \integral_{-\inf}^{\inf}f(t)exp(i*omega*t) dt
// Test that forward fft followed by inverse fft yields identity
TEST(FFTRVT_TEST, INVERSES)
{
    using spida::dcmplx;

    unsigned N = 32;
    spida::UniformGridRVT grid{N, -2, 2};
    unsigned nst = grid.getNst();

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0, 1.0);
    std::vector<double> in(N);
    std::vector<dcmplx> out(nst);
    std::vector<double> expect(N);
    for (unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::FFTRVT tr(grid);
    tr.T_To_ST(in, out);
    tr.ST_To_T(out, expect);

    EXPECT_LT(detail::relative_error(in, expect), 1e-10);
}

// FFT{exp(-(t/tp)^2)exp(-i*omega0*t}= tp*sqrt(pi)*exp(-tp^2*(omega-omega0)^2/4)
TEST(FFTRVT_TEST, GAUSST)
{
    using spida::dcmplx;
    using spida::PI;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;
    double minST = 1.10803e14;
    double maxST = 1.448963e16;

    spida::UniformGridRVT grid(nt, minT, maxT, minST, maxST);
    unsigned nst = grid.getNst();
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(nst);

    spida::FFTRVT transform(grid);
    spida::GaussT shape(grid, std::sqrt(I0), tp);
    shape.setFastPhase(omega0);
    auto y = shape.shapeRV();

    transform.T_To_ST(y, ysp);
    transform.ST_To_T(ysp, yinv);
    EXPECT_LT(detail::relative_error(y, yinv), 1e-6);

    std::vector<dcmplx> ysp_ex(grid.getNst(), 0.0);
    const std::vector<double>& omega = grid.getST();
    // y = f(t)*cos(i\omega0t) - > FFT{y} = (FFT{f(\omega - \omega0)}+FFT{f(\omega+\omega0)})/2
    // For real fields, fft taken over positive frequencies: FFT_real{y} =
    // FFT_real{f(\omega-\omega0)}/2
    for (size_t j = 0; j < grid.getNst(); j++)
        ysp_ex[j] = 0.5 * std::sqrt(I0) * tp * sqrt(PI) *
                    exp(-pow(tp, 2) * pow(omega[j] - omega0, 2) / 4.0);

    EXPECT_LT(detail::relative_error(ysp, ysp_ex), 1e-5);

    auto ycmplx = shape.shapeCV();
    // F{f(t)*exp(i*omega0*t)} = F{f(omega - omega0)}: expected complex-envelope spectrum
    std::vector<dcmplx> ysp_expectedCV(grid.getNst(), 0.0);
    for (size_t j = 0; j < grid.getNst(); j++)
        ysp_expectedCV[j] =
            std::sqrt(I0) * tp * sqrt(PI) * exp(-pow(tp, 2) * pow(omega[j] - omega0, 2) / 4.0);

    // Complex-valued transform output must match the expected Gaussian spectrum
    std::vector<dcmplx> ysp_outCV(grid.getNst(), 0.0);
    transform.CVT_To_ST(ycmplx, ysp_outCV);
    EXPECT_LT(detail::relative_error(ysp_expectedCV, ysp_outCV), 1e-5);
}

TEST(FFTRVT_TEST, COS)
{
    // F{cos(at)} = PI*(delta(omega-a) + delta(omega+a))
    unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    std::vector<double> in(N);
    spida::UniformGridRVT grid(N, 0.0, 2.0 * PI);
    unsigned nst = grid.getNst();
    std::vector<dcmplx> out(nst);

    const std::vector<double> t = grid.getT();
    for (size_t i = 0; i < t.size(); i++)
        in[i] = cos(8 * t[i]);

    spida::FFTRVT tr(grid);
    tr.T_To_ST(in, out);
    EXPECT_DOUBLE_EQ(out[8].real(), PI);
}

TEST(FFTRVT_TEST, DERIVATIVE_SIN)
{
    unsigned N = 32;

    using spida::dcmplx;
    using spida::PI;
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);

    // Need sin(tmin) = sin(tmax) for periodicity
    double tmin = 0.0;
    double tmax = 2.0 * PI;
    spida::UniformGridRVT grid(N, tmin, tmax);

    const std::vector<double> t = grid.getT();
    for (size_t i = 0; i < t.size(); i++)
        in[i] = sin(t[i]);
    for (size_t i = 0; i < t.size(); i++)
        expect[i] = cos(t[i]);

    spida::SpidaRVT spi(grid);
    spi.dT(in, out);
    EXPECT_LT(detail::relative_error(expect, out), 1e-6);

    std::vector<dcmplx> out1(grid.getNst());
    spida::FFTRVT tr(grid);
    tr.T_To_ST(in, out1);
    // sin(t) on [0, 2pi]: all spectral energy is in bin 1, purely imaginary
    EXPECT_NEAR(out1[1].real(), 0.0, 1e-10);
    EXPECT_NEAR(std::abs(out1[1]), spida::PI, 1e-10);
}

TEST(FFTCVT_TEST, ODD_N_THROWS)
{
    EXPECT_THROW(spida::FFTCVT(spida::UniformGridCVT{33u, -1.0, 1.0}), std::invalid_argument);
}

TEST(FFTRVT_TEST, ODD_N_THROWS)
{
    EXPECT_THROW(spida::FFTRVT(spida::UniformGridRVT{33u, -1.0, 1.0}), std::invalid_argument);
}

// ============================================================
//  P2: Non-power-of-2 even sizes (N=48) — never exercised
// ============================================================

TEST(FFTCVT_TEST, NON_POWER_OF_2_ROUND_TRIP)
{
    unsigned N = 48;
    using spida::dcmplx;
    std::default_random_engine gen(42);
    std::normal_distribution<double> dist;
    std::vector<dcmplx> in(N), sp(N), out(N);
    for (auto& v : in)
        v = dcmplx(dist(gen), dist(gen));

    spida::FFTCVT tr(spida::UniformGridCVT{N, -3.0, 3.0});
    tr.T_To_ST(in, sp);
    tr.ST_To_T(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

TEST(FFTRVT_TEST, NON_POWER_OF_2_ROUND_TRIP)
{
    unsigned N = 48;
    using spida::dcmplx;
    spida::UniformGridRVT grid(N, -3.0, 3.0);
    std::default_random_engine gen(43);
    std::normal_distribution<double> dist;
    std::vector<double> in(N), out(N);
    std::vector<dcmplx> sp(grid.getNst());
    for (auto& v : in)
        v = dist(gen);

    spida::FFTRVT tr(grid);
    tr.T_To_ST(in, sp);
    tr.ST_To_T(sp, out);
    EXPECT_LT(detail::relative_error(in, out), 1e-10);
}

// ============================================================
//  P2: FFTRVT band-limited round-trip with minI > 0
//  (fftRVT.cpp:85,92-93 — zero-fill loop for sub-band)
// ============================================================

TEST(FFTRVT_TEST, BAND_LIMITED_ROUND_TRIP_WITH_MIN_FREQ)
{
    // Choose minST > 0 so minI > 0: exercises the sub-band zero-fill loop
    // (fftRVT.cpp:85 and trailing zero-fill at lines 92-93).
    // Build the signal from the spectral domain to guarantee band-limitedness,
    // then verify ST -> T -> ST is the identity.
    unsigned N = 64;
    double minT = -4.0, maxT = 4.0;
    double dst = 2.0 * spida::PI / (maxT - minT);
    double minST = 3.0 * dst;
    double maxST = 10.0 * dst;
    spida::UniformGridRVT grid(N, minT, maxT, minST, maxST);
    EXPECT_GT(grid.getMinI(), 0u);

    using spida::dcmplx;
    std::default_random_engine gen(99);
    std::normal_distribution<double> dist;
    std::vector<dcmplx> sp_in(grid.getNst());
    for (auto& v : sp_in)
        v = dcmplx(dist(gen), dist(gen));

    // ST -> T uses the zero-fill; T -> ST should recover the original spectrum.
    std::vector<double> t_domain(N);
    std::vector<dcmplx> sp_out(grid.getNst());
    spida::FFTRVT tr(grid);
    tr.ST_To_T(sp_in, t_domain);
    tr.T_To_ST(t_domain, sp_out);
    EXPECT_LT(detail::relative_error(sp_in, sp_out), 1e-10);
}