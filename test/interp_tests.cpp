
#include <gtest/gtest.h>
#include <spida/helper/interp.h>
#include <pwutils/pwexcept.h>
#include <cmath>
#include <vector>

// --- tridisolve tests ---

TEST(TRIDISOLVE_TEST, KNOWN_SYSTEM)
{
    // Solve the 4x4 tridiagonal system:
    //  [ 2 -1  0  0 ] [x0]   [1]
    //  [-1  2 -1  0 ] [x1] = [0]
    //  [ 0 -1  2 -1 ] [x2]   [0]
    //  [ 0  0 -1  2 ] [x3]   [1]
    // Solution: x = [1, 1, 1, 1]
    std::vector<double> a = {-1.0, -1.0, -1.0};
    std::vector<double> b = { 2.0,  2.0,  2.0,  2.0};
    std::vector<double> c = {-1.0, -1.0, -1.0};
    std::vector<double> d = { 1.0,  0.0,  0.0,  1.0};
    std::vector<double> x;
    spida::tridisolve(a, b, c, d, x);
    ASSERT_EQ(x.size(), 4u);
    for (auto v : x)
        EXPECT_NEAR(v, 1.0, 1e-12);
}

TEST(TRIDISOLVE_TEST, SIZE_ERROR_SUBDIAG)
{
    std::vector<double> a = {1.0};              // size 1, too small
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), pw::Exception);
}

TEST(TRIDISOLVE_TEST, SIZE_ERROR_DIAG)
{
    std::vector<double> a = {1.0, 1.0};
    std::vector<double> b = {2.0, 2.0};         // size 2, too small
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), pw::Exception);
}

// --- LinearInterp tests ---

TEST(LINEAR_INTERP_TEST, EXACT_ON_LINEAR_FUNCTION)
{
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 2.0, 3.0};
    spida::LinearInterp interp(x, y);

    std::vector<double> xi = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    std::vector<double> yi = interp.eval(xi);
    ASSERT_EQ(yi.size(), xi.size());
    for (size_t i = 0; i < xi.size(); i++)
        EXPECT_NEAR(yi[i], xi[i], 1e-14);
}

TEST(LINEAR_INTERP_TEST, SINGLE_POINT_EVAL)
{
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 2.0, 4.0, 6.0};
    spida::LinearInterp interp(x, y);

    EXPECT_NEAR(interp.eval(0.0), 0.0, 1e-14);
    EXPECT_NEAR(interp.eval(1.5), 3.0, 1e-14);
    EXPECT_NEAR(interp.eval(3.0), 6.0, 1e-14);
}

TEST(LINEAR_INTERP_TEST, IN_PLACE_VECTOR_EVAL)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {1.0, 3.0, 5.0};
    spida::LinearInterp interp(x, y);

    std::vector<double> xi = {0.5, 1.0};
    std::vector<double> yi;
    interp.eval(xi, yi);
    ASSERT_EQ(yi.size(), 2u);
    EXPECT_NEAR(yi[0], 2.0, 1e-14);
    EXPECT_NEAR(yi[1], 3.0, 1e-14);
}

TEST(LINEAR_INTERP_TEST, ERROR_EMPTY_XINTERP)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    std::vector<double> xi;
    EXPECT_THROW(interp.eval(xi), pw::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_XINTERP_BELOW_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    EXPECT_THROW((void)interp.eval(-0.1), pw::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_XINTERP_ABOVE_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    EXPECT_THROW((void)interp.eval(2.1), pw::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_DATA_TOO_SMALL)
{
    std::vector<double> x = {1.0};
    std::vector<double> y = {1.0};
    EXPECT_THROW(spida::LinearInterp interp(x, y), pw::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_SIZE_MISMATCH)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0};
    EXPECT_THROW(spida::LinearInterp interp(x, y), pw::Exception);
}

// --- SplineInterp tests ---

TEST(SPLINE_INTERP_TEST, EXACT_ON_CUBIC_POLYNOMIAL)
{
    // Cubic spline interpolates cubic polynomials exactly.
    // Use f(x) = x^3 on [0, 3].
    unsigned N = 7;
    std::vector<double> x(N), y(N);
    for (unsigned i = 0; i < N; i++) {
        x[i] = i * 0.5;
        y[i] = x[i] * x[i] * x[i];
    }
    spida::SplineInterp interp(x, y);

    std::vector<double> xi = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75};
    std::vector<double> yi = interp.eval(xi);
    for (size_t i = 0; i < xi.size(); i++)
        EXPECT_NEAR(yi[i], xi[i] * xi[i] * xi[i], 1e-10);
}

TEST(SPLINE_INTERP_TEST, SMOOTH_FUNCTION_ACCURACY)
{
    // Spline on f(x) = sin(x) should be accurate.
    unsigned N = 20;
    double minx = 0.0;
    double maxx = 3.14159265358979;
    std::vector<double> x(N), y(N);
    for (unsigned i = 0; i < N; i++) {
        x[i] = minx + i * (maxx - minx) / (N - 1);
        y[i] = std::sin(x[i]);
    }
    spida::SplineInterp interp(x, y);

    std::vector<double> xi = {0.3, 0.7, 1.1, 1.5, 1.9, 2.3, 2.7};
    std::vector<double> yi = interp.eval(xi);
    for (size_t i = 0; i < xi.size(); i++)
        EXPECT_NEAR(yi[i], std::sin(xi[i]), 1e-5);
}

TEST(SPLINE_INTERP_TEST, SINGLE_POINT_EVAL)
{
    unsigned N = 10;
    std::vector<double> x(N), y(N);
    for (unsigned i = 0; i < N; i++) {
        x[i] = static_cast<double>(i);
        y[i] = x[i] * x[i];
    }
    spida::SplineInterp interp(x, y);
    EXPECT_NEAR(interp.eval(0.0), 0.0, 1e-10);
    EXPECT_NEAR(interp.eval(9.0), 81.0, 1e-10);
    EXPECT_NEAR(interp.eval(4.5), 4.5 * 4.5, 1e-3);
}

TEST(SPLINE_INTERP_TEST, IN_PLACE_VECTOR_EVAL)
{
    unsigned N = 6;
    std::vector<double> x(N), y(N);
    for (unsigned i = 0; i < N; i++) {
        x[i] = static_cast<double>(i);
        y[i] = static_cast<double>(i);
    }
    spida::SplineInterp interp(x, y);
    std::vector<double> xi = {1.0, 2.5, 4.0};
    std::vector<double> yi;
    interp.eval(xi, yi);
    ASSERT_EQ(yi.size(), 3u);
    EXPECT_NEAR(yi[0], 1.0, 1e-12);
    EXPECT_NEAR(yi[1], 2.5, 1e-8);
    EXPECT_NEAR(yi[2], 4.0, 1e-12);
}

TEST(SPLINE_INTERP_TEST, ERROR_XINTERP_OUT_OF_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 2.0, 3.0};
    spida::SplineInterp interp(x, y);
    EXPECT_THROW((void)interp.eval(-0.5), pw::Exception);
    EXPECT_THROW((void)interp.eval(3.5), pw::Exception);
}
