
#include <algorithm>
#include <cmath>
#include <iostream>

#include <gtest/gtest.h>
#include <spida/grid/besselR.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>
#include <spida/shape/shapeR.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeX.h>

TEST(SHAPE_TEST, GAUSST_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    spida::GaussT shape(grid, A, width);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
}

TEST(SHAPE_TEST, SECHT_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    spida::SechT shape(grid, A, width);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
}

TEST(SHAPE_TEST, AIRYT_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    double apod = 1.0;
    spida::AiryT shape(grid, A, width, apod);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
    EXPECT_DOUBLE_EQ(apod, shape.apodization());
}

TEST(SHAPE_TEST, SUPERGAUSST_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    double M = 2.0;
    spida::SuperGaussT shape(grid, A, width, M);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
    EXPECT_DOUBLE_EQ(M, shape.M());
}

TEST(SHAPE_TEST, BESSELT_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    double apod = 1.0;
    spida::BesselT shape(grid, A, width, apod);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
    EXPECT_DOUBLE_EQ(apod, shape.apodization());
}

TEST(SHAPE_TEST, GAUSST_COMPUTE)
{
    int N = 4096;
    double tmin = -240e-15;
    double tmax = 240e-15;
    double freq_min = 1.10803e14;
    double freq_max = 1.448963e16;
    spida::UniformGridRVT grid(N, tmin, tmax, freq_min, freq_max);

    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;

    spida::GaussT shape(grid, std::sqrt(I0), tp);
    shape.setFastPhase(car_freq);
    auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y), std::end(y));
    auto itmin = std::min_element(std::begin(y), std::end(y));
    EXPECT_NEAR(223606797.74997896, *itmax, 1e-6);
    EXPECT_NEAR(-200519096.49044815, *itmin, 1e-6);
    EXPECT_NEAR(90533935.926272079, y[(N - N / 16) / 2], 1e-6);
    EXPECT_NEAR(90533935.926272079, y[(N + N / 16) / 2], 1e-6);
}

TEST(SHAPE_TEST, SECHT_COMPUTE)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;

    spida::UniformGridRVT grid(nt, -240e-15, 240e-15, 1.10803e14, 1.448963e16);
    spida::SechT shape(grid, std::sqrt(I0), tp);
    shape.setFastPhase(car_freq);
    auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y), std::end(y));
    auto itmin = std::min_element(std::begin(y), std::end(y));

    EXPECT_NEAR(223606797.74997896, *itmax, 1e-6);
    EXPECT_NEAR(-211808302.73070601, *itmin, 1e-6);
    EXPECT_NEAR(122726544.58500151, y[(nt - nt / 16) / 2], 1e-6);
    EXPECT_NEAR(122726544.58500151, y[(nt + nt / 16) / 2], 1e-6);
}

TEST(SHAPE_TEST, AIRYT_COMPUTE)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -320e-15;
    double maxT = 120e-15;
    double apod = 0.2;

    spida::UniformGridRVT grid(nt, minT, maxT, 1.10803e14, 1.448963e16);
    spida::AiryT shape(grid, std::sqrt(I0), tp, apod);
    shape.setFastPhase(omega0);
    auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y), std::end(y));
    auto itmin = std::min_element(std::begin(y), std::end(y));

    EXPECT_NEAR(111197223.83176377, *itmax, 1e-6);
    EXPECT_NEAR(-115045662.1447591, *itmin, 1e-6);
    EXPECT_NEAR(8518347.9235188272, y[(nt - nt / 16) / 2], 1e-6);
    EXPECT_NEAR(-18030265.656660724, y[(nt + nt / 16) / 2], 1e-6);
}

TEST(SHAPE_TEST, BESSELT_COMPUTE)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;
    double apod = 0.5;

    spida::UniformGridRVT grid(nt, minT, maxT, 1.10803e14, 1.448963e16);
    spida::BesselT shape(grid, std::sqrt(I0), tp, apod);
    shape.setFastPhase(omega0);
    auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y), std::end(y));
    auto itmin = std::min_element(std::begin(y), std::end(y));

    EXPECT_NEAR(223606797.74997896, *itmax, 1e-6);
    EXPECT_NEAR(-185314259.84787187, *itmin, 1e-6);
    EXPECT_NEAR(46643813.615054831, y[(nt - nt / 16) / 2], 1e-6);
    EXPECT_NEAR(46643813.615054831, y[(nt + nt / 16) / 2], 1e-6);
}

// ============================================================
//  Additional ShapeT tests — setters, envelope, shapeCV, chirp
// ============================================================

TEST(SHAPE_TEST, GAUSST_SETTERS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    spida::GaussT shape(grid, 1.0, 1.0);

    shape.setAmplitude(3.0);
    shape.setWidth(2.0);
    shape.setChirp(0.5);
    shape.setOffset(0.1);
    shape.setSlowPhase(0.2);
    shape.setFastPhase(0.3);

    EXPECT_DOUBLE_EQ(shape.amplitude(), 3.0);
    EXPECT_DOUBLE_EQ(shape.width(), 2.0);
    EXPECT_DOUBLE_EQ(shape.chirp(), 0.5);
    EXPECT_DOUBLE_EQ(shape.offset(), 0.1);
    EXPECT_DOUBLE_EQ(shape.slowPhase(), 0.2);
    EXPECT_DOUBLE_EQ(shape.fastPhase(), 0.3);
}

TEST(SHAPE_TEST, GAUSST_SHAPECV_SIZE)
{
    unsigned nt = 64;
    spida::UniformGridRVT grid(nt, -5.0, 5.0);
    spida::GaussT shape(grid, 1.0, 1.0);
    EXPECT_EQ(shape.shapeCV().size(), nt);
    EXPECT_EQ(shape.getT().size(), nt);
}

TEST(SHAPE_TEST, GAUSST_ENVELOPE_EQUALS_SHAPECV_WHEN_FAST_PHASE_ZERO)
{
    // With omega0=0 the fast phase factor is exp(0)=1, so shapeCV == envelope
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    spida::GaussT shape(grid, 2.0, 1.0); // omega0 defaults to 0

    auto env = shape.envelope();
    auto cv = shape.shapeCV();
    ASSERT_EQ(env.size(), cv.size());
    for (size_t i = 0; i < env.size(); i++) {
        EXPECT_NEAR(env[i].real(), cv[i].real(), 1e-14);
        EXPECT_NEAR(env[i].imag(), cv[i].imag(), 1e-14);
    }
}

TEST(SHAPE_TEST, GAUSST_ENVELOPE_PEAK_AT_CENTER)
{
    // Grid with 64 points from -5 to 5; dx = 10/64, so t[32] = -5 + 32*(10/64) = 0
    unsigned nt = 64;
    spida::UniformGridRVT grid(nt, -5.0, 5.0);
    double A = 3.0;
    spida::GaussT shape(grid, A, 1.0);

    auto env = shape.envelope();
    // envelope peak should equal A (at t=0, compute(0)=exp(0)=1, slowPhase=0)
    EXPECT_NEAR(env[nt / 2].real(), A, 1e-12);
    EXPECT_NEAR(env[nt / 2].imag(), 0.0, 1e-12);
}

TEST(SHAPE_TEST, GAUSST_CHIRP_MAKES_ENVELOPE_COMPLEX)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    spida::GaussT shape(grid, 1.0, 1.0);
    shape.setChirp(1.0);

    auto env = shape.envelope();
    // At t != 0, slowPhaseFactor = exp(-ii*chirp*t^2) has non-zero imaginary part
    double max_imag = 0.0;
    for (const auto& v : env)
        max_imag = std::max(max_imag, std::abs(v.imag()));
    EXPECT_GT(max_imag, 0.0);
}

TEST(SHAPE_TEST, SECHT_SHAPECV_PEAK)
{
    unsigned nt = 64;
    spida::UniformGridRVT grid(nt, -5.0, 5.0);
    double A = 2.0;
    spida::SechT shape(grid, A, 1.0); // compute(0) = 1/cosh(0) = 1

    auto cv = shape.shapeCV();
    EXPECT_NEAR(cv[nt / 2].real(), A, 1e-12);
}

TEST(SHAPE_TEST, SUPERGAUSST_COMPUTE)
{
    unsigned nt = 64;
    spida::UniformGridRVT grid(nt, -5.0, 5.0);
    double A = 1.0, tp = 1.0, M = 2.0;
    spida::SuperGaussT shape(grid, A, tp, M);

    auto rv = shape.shapeRV();
    // Peak at t=0: compute(0) = exp(0) = 1
    EXPECT_NEAR(rv[nt / 2], A, 1e-12);
    // Values are positive (super-Gaussian is non-negative)
    for (double v : rv)
        EXPECT_GE(v, 0.0);
}

// ============================================================
//  ShapeX tests — ExpX, CosX, SinX
// ============================================================
//  Grid: UniformGridRVX(4, 0.0, 4.0)
//    dx = (4-0)/4 = 1 -> x = {0, 1, 2, 3}

TEST(SHAPEX_TEST, EXPX_SET_VARIABLES)
{
    spida::UniformGridRVX grid(4, 0.0, 4.0);
    spida::ExpX shape(grid, 2.0, 0.5);

    EXPECT_DOUBLE_EQ(shape.amplitude(), 2.0);
    EXPECT_DOUBLE_EQ(shape.width(), 0.5);
    EXPECT_DOUBLE_EQ(shape.offset(), 0.0);
    EXPECT_DOUBLE_EQ(shape.phase(), 0.0);
    EXPECT_EQ(shape.getX().size(), 4u);
}

TEST(SHAPEX_TEST, EXPX_SHAPERV_VALUES)
{
    // shapeRV()[i] = A * exp((x[i]-offset)/w0) * cos(phi0)
    // A=1, w0=1, offset=0, phi0=0 -> exp(x[i])
    spida::UniformGridRVX grid(4, 0.0, 4.0);
    spida::ExpX shape(grid, 1.0, 1.0);
    auto rv = shape.shapeRV();

    ASSERT_EQ(rv.size(), 4u);
    EXPECT_NEAR(rv[0], std::exp(0.0), 1e-12);
    EXPECT_NEAR(rv[1], std::exp(1.0), 1e-12);
    EXPECT_NEAR(rv[2], std::exp(2.0), 1e-12);
    EXPECT_NEAR(rv[3], std::exp(3.0), 1e-12);
}

TEST(SHAPEX_TEST, COSX_SHAPERV_VALUES)
{
    // shapeRV()[i] = A * cos(x[i]/w0) with A=1, w0=1, x={0,1,2,3}
    spida::UniformGridRVX grid(4, 0.0, 4.0);
    spida::CosX shape(grid, 1.0, 1.0);
    auto rv = shape.shapeRV();

    ASSERT_EQ(rv.size(), 4u);
    EXPECT_NEAR(rv[0], std::cos(0.0), 1e-12);
    EXPECT_NEAR(rv[1], std::cos(1.0), 1e-12);
    EXPECT_NEAR(rv[2], std::cos(2.0), 1e-12);
    EXPECT_NEAR(rv[3], std::cos(3.0), 1e-12);
}

TEST(SHAPEX_TEST, SINX_SHAPERV_VALUES)
{
    // shapeRV()[i] = A * sin(x[i]/w0) with A=1, w0=1, x={0,1,2,3}
    spida::UniformGridRVX grid(4, 0.0, 4.0);
    spida::SinX shape(grid, 1.0, 1.0);
    auto rv = shape.shapeRV();

    ASSERT_EQ(rv.size(), 4u);
    EXPECT_NEAR(rv[0], std::sin(0.0), 1e-12);
    EXPECT_NEAR(rv[1], std::sin(1.0), 1e-12);
    EXPECT_NEAR(rv[2], std::sin(2.0), 1e-12);
    EXPECT_NEAR(rv[3], std::sin(3.0), 1e-12);
}

TEST(SHAPEX_TEST, SHAPERV_MATCHES_REAL_PART_OF_SHAPECV)
{
    // With zero phase, real(shapeCV()[i]) == shapeRV()[i]
    spida::UniformGridRVX grid(8, -4.0, 4.0);
    spida::CosX shape(grid, 1.0, 1.0);
    auto rv = shape.shapeRV();
    auto cv = shape.shapeCV();

    ASSERT_EQ(rv.size(), cv.size());
    for (size_t i = 0; i < rv.size(); i++)
        EXPECT_NEAR(cv[i].real(), rv[i], 1e-12);
}

TEST(SHAPEX_TEST, SHAPECV_IMAG_ZERO_WITHOUT_PHASE)
{
    spida::UniformGridRVX grid(8, -4.0, 4.0);
    spida::SinX shape(grid, 1.0, 1.0);
    auto cv = shape.shapeCV();
    for (const auto& v : cv)
        EXPECT_NEAR(v.imag(), 0.0, 1e-12);
}

TEST(SHAPEX_TEST, PHASE_ROTATES_SHAPECV)
{
    // With phi0=pi/2: shapeCV()[i] = cos(x[i]) * exp(ii*pi/2) = ii*cos(x[i])
    spida::UniformGridRVX grid(4, 0.0, 4.0);
    spida::CosX shape(grid, 1.0, 1.0);
    shape.setPhase(spida::PI / 2.0);
    auto cv = shape.shapeCV();

    // real = cos(x)*cos(pi/2) = 0, imag = cos(x)*sin(pi/2) = cos(x)
    EXPECT_NEAR(cv[0].real(), 0.0, 1e-12);
    EXPECT_NEAR(cv[0].imag(), std::cos(0.0), 1e-12);
    EXPECT_NEAR(cv[1].imag(), std::cos(1.0), 1e-12);
}

TEST(SHAPEX_TEST, OFFSET_SHIFTS_ARGUMENT)
{
    // CosX at x=1 with offset=1 equals CosX at x=0 without offset
    spida::UniformGridRVX grid(4, 0.0, 4.0); // x = {0, 1, 2, 3}

    spida::CosX no_offset(grid, 1.0, 1.0);
    spida::CosX with_offset(grid, 1.0, 1.0);
    with_offset.setOffset(1.0);

    auto rv_none = no_offset.shapeRV();  // cos(x[i])
    auto rv_off = with_offset.shapeRV(); // cos(x[i] - 1)

    // rv_off[1] = cos((1-1)/1) = cos(0) = 1 == rv_none[0]
    EXPECT_NEAR(rv_off[1], rv_none[0], 1e-12);
    EXPECT_NEAR(rv_off[2], rv_none[1], 1e-12);
}

TEST(SHAPEX_TEST, SETTERS)
{
    spida::UniformGridRVX grid(4, 0.0, 4.0);
    spida::SinX shape(grid, 1.0, 1.0);

    shape.setAmplitude(5.0);
    shape.setWidth(2.0);
    shape.setOffset(0.5);
    shape.setPhase(0.3);

    EXPECT_DOUBLE_EQ(shape.amplitude(), 5.0);
    EXPECT_DOUBLE_EQ(shape.width(), 2.0);
    EXPECT_DOUBLE_EQ(shape.offset(), 0.5);
    EXPECT_DOUBLE_EQ(shape.phase(), 0.3);
}

// ============================================================
//  ShapeR tests — GaussR
// ============================================================

TEST(SHAPER_TEST, GAUSSR_SET_VARIABLES)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 2.0, 1.0);

    EXPECT_DOUBLE_EQ(shape.amplitude(), 2.0);
    EXPECT_DOUBLE_EQ(shape.width(), 1.0);
    EXPECT_EQ(shape.getR().size(), 16u);
}

TEST(SHAPER_TEST, GAUSSR_SHAPERV_SIZE_AND_POSITIVE)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 2.0);
    auto rv = shape.shapeRV();

    EXPECT_EQ(rv.size(), 16u);
    for (double v : rv)
        EXPECT_GT(v, 0.0); // GaussR is always positive
}

TEST(SHAPER_TEST, GAUSSR_SHAPERV_DECREASING)
{
    // Bessel grid r points are increasing; Gaussian decreases with r
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 2.0);
    auto rv = shape.shapeRV();

    EXPECT_TRUE(std::is_sorted(rv.begin(), rv.end(), std::greater<double>{}));
}

TEST(SHAPER_TEST, GAUSSR_SHAPECV_MATCHES_SHAPERV_WITHOUT_FOCUS)
{
    // Without focus, shapeCV is real-valued (imaginary ≈ 0)
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 2.0);
    auto rv = shape.shapeRV();
    auto cv = shape.shapeCV();

    ASSERT_EQ(rv.size(), cv.size());
    for (size_t i = 0; i < rv.size(); i++) {
        EXPECT_NEAR(cv[i].real(), rv[i], 1e-12);
        EXPECT_NEAR(cv[i].imag(), 0.0, 1e-12);
    }
}

TEST(SHAPER_TEST, GAUSSR_SET_FOCUS_ADDS_PHASE)
{
    // With focus f, shapeCV()[i] = shapeRV()[i] * exp(-ii*r[i]^2/f)
    // -> imaginary part is non-zero for r > 0
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 2.0);
    shape.setFocus(10.0);

    EXPECT_NEAR(shape.focus(), 10.0, 1e-15);

    auto cv = shape.shapeCV();
    double max_imag = 0.0;
    for (const auto& v : cv)
        max_imag = std::max(max_imag, std::abs(v.imag()));
    EXPECT_GT(max_imag, 0.0);
}

TEST(SHAPER_TEST, GAUSSR_SET_FOCUS_ZERO_THROWS)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 2.0);
    EXPECT_THROW(shape.setFocus(0.0), std::domain_error);
}

TEST(SHAPER_TEST, GAUSSR_SET_FOCUS_NEGATIVE_THROWS)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 2.0);
    EXPECT_THROW(shape.setFocus(-5.0), std::domain_error);
}

TEST(SHAPER_TEST, GAUSSR_SETTERS)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 1.0);

    shape.setAmplitude(4.0);
    shape.setWidth(3.0);

    EXPECT_DOUBLE_EQ(shape.amplitude(), 4.0);
    EXPECT_DOUBLE_EQ(shape.width(), 3.0);
}

TEST(SHAPER_TEST, GAUSSR_SHAPERV_SCALES_WITH_AMPLITUDE)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape1(grid, 1.0, 2.0);
    spida::GaussR shape2(grid, 3.0, 2.0);

    auto rv1 = shape1.shapeRV();
    auto rv2 = shape2.shapeRV();

    for (size_t i = 0; i < rv1.size(); i++)
        EXPECT_NEAR(rv2[i], 3.0 * rv1[i], 1e-12);
}

// ============================================================
//  Zero-width constructor guard tests — ShapeT
// ============================================================

TEST(SHAPE_TEST, GAUSST_ZERO_WIDTH_THROWS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    EXPECT_THROW(spida::GaussT(grid, 1.0, 0.0), std::domain_error);
}

TEST(SHAPE_TEST, SECHT_ZERO_WIDTH_THROWS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    EXPECT_THROW(spida::SechT(grid, 1.0, 0.0), std::domain_error);
}

TEST(SHAPE_TEST, AIRYT_ZERO_WIDTH_THROWS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    EXPECT_THROW(spida::AiryT(grid, 1.0, 0.0, 0.1), std::domain_error);
}

TEST(SHAPE_TEST, SUPERGAUSST_ZERO_WIDTH_THROWS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    EXPECT_THROW(spida::SuperGaussT(grid, 1.0, 0.0, 2.0), std::domain_error);
}

TEST(SHAPE_TEST, BESSELT_ZERO_WIDTH_THROWS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    EXPECT_THROW(spida::BesselT(grid, 1.0, 0.0, 0.5), std::domain_error);
}

// ============================================================
//  setWidth zero guard — ShapeT
// ============================================================

TEST(SHAPE_TEST, GAUSST_SETWIDTH_ZERO_THROWS)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    spida::GaussT shape(grid, 1.0, 1.0);
    EXPECT_THROW(shape.setWidth(0.0), std::domain_error);
}

TEST(SHAPE_TEST, GAUSST_SETWIDTH_VALID_UPDATES)
{
    spida::UniformGridRVT grid(64, -5.0, 5.0);
    spida::GaussT shape(grid, 1.0, 1.0);
    shape.setWidth(2.0);
    EXPECT_DOUBLE_EQ(shape.width(), 2.0);
}

// ============================================================
//  Zero-width constructor guard tests — ShapeR
// ============================================================

TEST(SHAPER_TEST, GAUSSR_ZERO_WIDTH_THROWS)
{
    spida::BesselRootGridR grid(16, 5.0);
    EXPECT_THROW(spida::GaussR(grid, 1.0, 0.0), std::domain_error);
}

// ============================================================
//  setWidth zero guard — ShapeR
// ============================================================

TEST(SHAPER_TEST, GAUSSR_SETWIDTH_ZERO_THROWS)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 1.0);
    EXPECT_THROW(shape.setWidth(0.0), std::domain_error);
}

TEST(SHAPER_TEST, GAUSSR_SETWIDTH_VALID_UPDATES)
{
    spida::BesselRootGridR grid(16, 5.0);
    spida::GaussR shape(grid, 1.0, 1.0);
    shape.setWidth(2.0);
    EXPECT_DOUBLE_EQ(shape.width(), 2.0);
}

// ============================================================
//  Zero-width constructor guard tests — ShapeX
// ============================================================

TEST(SHAPEX_TEST, SHAPEX_ZERO_WIDTH_THROWS)
{
    // SinX exercises the ShapeX base-class constructor guard shared by all
    // ShapeX concrete types (ExpX, SinX, CosX).
    spida::UniformGridRVX grid(8, -4.0, 4.0);
    EXPECT_THROW(spida::SinX(grid, 1.0, 0.0), std::domain_error);
}

// ============================================================
//  setWidth zero guard — ShapeX
// ============================================================

TEST(SHAPEX_TEST, SHAPEX_SETWIDTH_ZERO_THROWS)
{
    spida::UniformGridRVX grid(8, -4.0, 4.0);
    spida::SinX shape(grid, 1.0, 1.0);
    EXPECT_THROW(shape.setWidth(0.0), std::domain_error);
}

TEST(SHAPEX_TEST, SHAPEX_SETWIDTH_VALID_UPDATES)
{
    spida::UniformGridRVX grid(8, -4.0, 4.0);
    spida::SinX shape(grid, 1.0, 1.0);
    shape.setWidth(2.0);
    EXPECT_DOUBLE_EQ(shape.width(), 2.0);
}