
#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <spida/grid/chebX.h>
#include <spida/grid/uniformCVT.h>
#include <spida/grid/uniformCVX.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>

// --- UniformGridCVX tests ---

TEST(UNIFORM_GRID_CVX_TEST, BASIC_PROPERTIES)
{
    unsigned N = 32;
    double minX = -3.0;
    double maxX = 3.0;
    spida::UniformGridCVX grid(N, minX, maxX);

    EXPECT_EQ(grid.getNx(), N);
    EXPECT_DOUBLE_EQ(grid.getMinX(), minX);
    EXPECT_DOUBLE_EQ(grid.getMaxX(), maxX);
    EXPECT_DOUBLE_EQ(grid.getLX(), maxX - minX);
    EXPECT_EQ(grid.getX().size(), N);
    EXPECT_EQ(grid.getSX().size(), N);
}

TEST(UNIFORM_GRID_CVX_TEST, GRID_SPACING)
{
    unsigned N = 8;
    double minX = 0.0;
    double maxX = 4.0;
    spida::UniformGridCVX grid(N, minX, maxX);
    const auto& x = grid.getX();
    double expected_dx = (maxX - minX) / N;
    for (unsigned i = 1; i < N; i++)
        EXPECT_NEAR(x[i] - x[i - 1], expected_dx, 1e-14);
}

TEST(UNIFORM_GRID_CVX_TEST, SPECTRAL_GRID_RANGE)
{
    unsigned N = 16;
    double minX = -spida::PI;
    double maxX = spida::PI;
    spida::UniformGridCVX grid(N, minX, maxX);

    EXPECT_NEAR(grid.getMinSX(), -(N * spida::PI) / (maxX - minX), 1e-12);
    EXPECT_NEAR(grid.getMaxSX(), (N * spida::PI) / (maxX - minX), 1e-12);
}

TEST(UNIFORM_GRID_CVX_TEST, FREQSHIFT_REAL_ROUNDTRIP)
{
    // freqshift applied twice does not restore original for odd-style shifts,
    // but the output size and first element can be validated.
    unsigned N = 8;
    spida::UniformGridCVX grid(N, -1.0, 1.0);
    const auto& sx = grid.getSX();
    auto shifted = grid.freqshift(sx);
    EXPECT_EQ(shifted.size(), N);
    // After freqshift the negative frequencies come first, verify monotonicity
    for (unsigned i = 1; i < shifted.size(); i++)
        EXPECT_LE(shifted[i - 1], shifted[i] + 1e-14);
}

TEST(UNIFORM_GRID_CVX_TEST, FREQSHIFT_INPLACE_REAL)
{
    unsigned N = 8;
    spida::UniformGridCVX grid(N, -1.0, 1.0);
    const auto& sx = grid.getSX();
    std::vector<double> out(N);
    grid.freqshift(sx, out);
    auto out2 = grid.freqshift(sx);
    EXPECT_EQ(out, out2);
}

TEST(UNIFORM_GRID_CVX_TEST, ERROR_INVALID_RANGE)
{
    EXPECT_THROW(spida::UniformGridCVX(16, 1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(spida::UniformGridCVX(16, 2.0, 2.0), std::invalid_argument);
}

// --- UniformGridRVX tests ---

TEST(UNIFORM_GRID_RVX_TEST, BASIC_PROPERTIES)
{
    unsigned N = 32;
    double minX = -spida::PI;
    double maxX = spida::PI;
    spida::UniformGridRVX grid(N, minX, maxX);

    EXPECT_EQ(grid.getNx(), N);
    EXPECT_DOUBLE_EQ(grid.getMinX(), minX);
    EXPECT_DOUBLE_EQ(grid.getMaxX(), maxX);
    EXPECT_EQ(grid.getX().size(), N);
    EXPECT_EQ(grid.getNsx(), N / 2 + 1);
    EXPECT_EQ(grid.getSX().size(), N / 2 + 1);
}

TEST(UNIFORM_GRID_RVX_TEST, SPECTRAL_GRID_NON_NEGATIVE)
{
    unsigned N = 16;
    spida::UniformGridRVX grid(N, -1.0, 1.0);
    const auto& sx = grid.getSX();
    for (double val : sx)
        EXPECT_GE(val, 0.0);
}

TEST(UNIFORM_GRID_RVX_TEST, SPECTRAL_MIN_MAX)
{
    unsigned N = 16;
    double minX = -2.0;
    double maxX = 2.0;
    spida::UniformGridRVX grid(N, minX, maxX);
    EXPECT_DOUBLE_EQ(grid.getMinSX(), 0.0);
    EXPECT_NEAR(grid.getMaxSX(), N * spida::PI / (maxX - minX), 1e-12);
}

// --- UniformGridCVT tests ---

TEST(UNIFORM_GRID_CVT_TEST, BASIC_PROPERTIES)
{
    unsigned N = 64;
    double minT = -5.0;
    double maxT = 5.0;
    spida::UniformGridCVT grid(N, minT, maxT);

    EXPECT_EQ(grid.getNt(), N);
    EXPECT_DOUBLE_EQ(grid.getMinT(), minT);
    EXPECT_DOUBLE_EQ(grid.getMaxT(), maxT);
    EXPECT_DOUBLE_EQ(grid.getLT(), maxT - minT);
    EXPECT_EQ(grid.getT().size(), N);
    EXPECT_EQ(grid.getST().size(), N);
    EXPECT_EQ(grid.getNst(), N);
}

TEST(UNIFORM_GRID_CVT_TEST, SPECTRAL_RANGE)
{
    unsigned N = 64;
    double minT = -1.0;
    double maxT = 1.0;
    spida::UniformGridCVT grid(N, minT, maxT);

    EXPECT_NEAR(grid.getMinST(), -(N * spida::PI) / (maxT - minT), 1e-12);
    EXPECT_NEAR(grid.getMaxST(), (N * spida::PI) / (maxT - minT), 1e-12);
}

TEST(UNIFORM_GRID_CVT_TEST, ERROR_ODD_NT)
{
    EXPECT_THROW(spida::UniformGridCVT(33, -1.0, 1.0), std::invalid_argument);
}

TEST(UNIFORM_GRID_CVT_TEST, FREQSHIFT_REAL_MONOTONE)
{
    unsigned N = 16;
    spida::UniformGridCVT grid(N, -1.0, 1.0);
    const auto& st = grid.getST();
    auto shifted = grid.freqshift(st);
    EXPECT_EQ(shifted.size(), N);
    for (unsigned i = 1; i < shifted.size(); i++)
        EXPECT_LE(shifted[i - 1], shifted[i] + 1e-14);
}

TEST(UNIFORM_GRID_CVT_TEST, FREQSHIFT_INPLACE_MATCHES_RETURN)
{
    unsigned N = 16;
    spida::UniformGridCVT grid(N, -1.0, 1.0);
    const auto& st = grid.getST();
    std::vector<double> out(N);
    grid.freqshift(st, out);
    auto out2 = grid.freqshift(st);
    EXPECT_EQ(out, out2);
}

// --- UniformGridRVT tests ---

TEST(UNIFORM_GRID_RVT_TEST, BASIC_PROPERTIES)
{
    unsigned N = 128;
    double minT = -240e-15;
    double maxT = 240e-15;
    spida::UniformGridRVT grid(N, minT, maxT);

    EXPECT_EQ(grid.getNt(), N);
    EXPECT_DOUBLE_EQ(grid.getMinT(), minT);
    EXPECT_DOUBLE_EQ(grid.getMaxT(), maxT);
    EXPECT_EQ(grid.getNst(), N / 2 + 1);
    EXPECT_EQ(grid.getST().size(), N / 2 + 1);
}

TEST(UNIFORM_GRID_RVT_TEST, SPECTRAL_NON_NEGATIVE)
{
    unsigned N = 64;
    spida::UniformGridRVT grid(N, -1.0, 1.0);
    const auto& st = grid.getST();
    for (double val : st)
        EXPECT_GE(val, 0.0);
}

TEST(UNIFORM_GRID_RVT_TEST, FREQ_RANGE_CONSTRUCTOR)
{
    unsigned N = 4096;
    double minT = -240e-15;
    double maxT = 240e-15;
    double minST = 1.1e14;
    double maxST = 1.4e16;
    spida::UniformGridRVT grid(N, minT, maxT, minST, maxST);

    // Grid snaps to nearest bin, so actual min/max are within one bin-width of requested
    double dst = 2.0 * spida::PI / (maxT - minT);
    EXPECT_NEAR(grid.getMinST(), minST, dst);
    EXPECT_NEAR(grid.getMaxST(), maxST, dst);
    EXPECT_LT(grid.getNst(), N / 2 + 1);
    EXPECT_GT(grid.getNst(), 0u);
}

TEST(UNIFORM_GRID_RVT_TEST, FREQ_TO_INDEX_AND_BACK)
{
    unsigned N = 256;
    double minT = -1.0;
    double maxT = 1.0;
    spida::UniformGridRVT grid(N, minT, maxT);

    unsigned idx = grid.freqToIndx(0.0);
    EXPECT_EQ(idx, 0u);

    double freq0 = grid.indxToFreq(0);
    EXPECT_NEAR(freq0, 0.0, 1e-14);

    double freqMax = grid.indxToFreq(N / 2);
    EXPECT_NEAR(freqMax, grid.maxPossibleFreq(), 1e-12);
}

TEST(UNIFORM_GRID_RVT_TEST, ERROR_INDX_OUT_OF_RANGE)
{
    unsigned N = 64;
    spida::UniformGridRVT grid(N, -1.0, 1.0);
    EXPECT_THROW((void) grid.indxToFreq(N / 2 + 1), std::domain_error);
}

TEST(UNIFORM_GRID_RVT_TEST, ERROR_FREQ_NEGATIVE)
{
    unsigned N = 64;
    spida::UniformGridRVT grid(N, -1.0, 1.0);
    EXPECT_THROW((void) grid.freqToIndx(-1.0), std::domain_error);
}

TEST(UNIFORM_GRID_RVT_TEST, ERROR_INVALID_FREQ_RANGE)
{
    // maxST > maxPossibleFreq
    unsigned N = 64;
    double minT = -1.0;
    double maxT = 1.0;
    double minST = 0.0;
    double maxST = 1e20;
    EXPECT_THROW(spida::UniformGridRVT grid(N, minT, maxT, minST, maxST), std::domain_error);
}

// --- ChebGridX tests ---

TEST(CHEB_GRID_X_TEST, EXTREMA_PROPERTIES)
{
    int N = 16;
    double minX = -1.0;
    double maxX = 1.0;
    spida::ChebExtremaGridX grid(N, minX, maxX);

    EXPECT_EQ(grid.getNx(), static_cast<unsigned>(N));
    EXPECT_DOUBLE_EQ(grid.getMinX(), minX);
    EXPECT_DOUBLE_EQ(grid.getMaxX(), maxX);
    EXPECT_EQ(grid.getX().size(), static_cast<size_t>(N));
    EXPECT_EQ(grid.getSX().size(), static_cast<size_t>(N));
    EXPECT_DOUBLE_EQ(grid.getMinSX(), 0.0);
    EXPECT_NEAR(grid.getMaxSX(), spida::PI, 1e-14);
}

TEST(CHEB_GRID_X_TEST, EXTREMA_ENDPOINTS)
{
    int N = 16;
    spida::ChebExtremaGridX grid(N, -1.0, 1.0);
    const auto& x = grid.getX();
    EXPECT_NEAR(std::abs(x.front()), 1.0, 1e-12);
    EXPECT_NEAR(std::abs(x.back()), 1.0, 1e-12);
}

TEST(CHEB_GRID_X_TEST, ROOT_PROPERTIES)
{
    int N = 16;
    spida::ChebRootGridX grid(N, -1.0, 1.0);

    EXPECT_EQ(grid.getNx(), static_cast<unsigned>(N));
    EXPECT_EQ(grid.getX().size(), static_cast<size_t>(N));
    EXPECT_EQ(grid.getSX().size(), static_cast<size_t>(N));
}

TEST(CHEB_GRID_X_TEST, ROOT_INTERIOR_POINTS)
{
    int N = 8;
    spida::ChebRootGridX grid(N, -1.0, 1.0);
    const auto& x = grid.getX();
    // Chebyshev root nodes are strictly inside (-1, 1)
    for (double xi : x) {
        EXPECT_GT(xi, -1.0);
        EXPECT_LT(xi, 1.0);
    }
}
