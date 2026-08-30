
#include <cmath>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>
#include <spida/helper/interp.h>
#include <utils/except.h>

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
    std::vector<double> b = {2.0, 2.0, 2.0, 2.0};
    std::vector<double> c = {-1.0, -1.0, -1.0};
    std::vector<double> d = {1.0, 0.0, 0.0, 1.0};
    std::vector<double> x;
    spida::tridisolve(a, b, c, d, x);
    ASSERT_EQ(x.size(), 4u);
    for (auto v : x)
        EXPECT_NEAR(v, 1.0, 1e-12);
}

TEST(TRIDISOLVE_TEST, SIZE_ERROR_SUBDIAG)
{
    std::vector<double> a = {1.0}; // size 1, too small
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
}

TEST(TRIDISOLVE_TEST, SIZE_ERROR_DIAG)
{
    std::vector<double> a = {1.0, 1.0};
    std::vector<double> b = {2.0, 2.0}; // size 2, too small
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
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
    EXPECT_THROW((void) interp.eval(xi), detail::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_XINTERP_BELOW_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    EXPECT_THROW((void) interp.eval(-0.1), detail::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_XINTERP_ABOVE_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    EXPECT_THROW((void) interp.eval(2.1), detail::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_DATA_TOO_SMALL)
{
    std::vector<double> x = {1.0};
    std::vector<double> y = {1.0};
    EXPECT_THROW(spida::LinearInterp interp(x, y), detail::Exception);
}

TEST(LINEAR_INTERP_TEST, ERROR_SIZE_MISMATCH)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0};
    EXPECT_THROW(spida::LinearInterp interp(x, y), detail::Exception);
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
    double maxx = std::numbers::pi;
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
    EXPECT_NEAR(interp.eval(4.5), 4.5 * 4.5, 1e-10);
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

TEST(SPLINE_INTERP_TEST, ERROR_XINTERP_BELOW_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 2.0, 3.0};
    spida::SplineInterp interp(x, y);
    EXPECT_THROW((void) interp.eval(-0.5), detail::Exception);
}

TEST(SPLINE_INTERP_TEST, ERROR_XINTERP_ABOVE_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 2.0, 3.0};
    spida::SplineInterp interp(x, y);
    EXPECT_THROW((void) interp.eval(3.5), detail::Exception);
}

// --- tridisolve guard-clause tests (P1 gap coverage) ---

// c.size() < 2: superdiagonal vector has only one element (line 110-112)
TEST(TRIDISOLVE_TEST, SIZE_ERROR_SUPERDIAG)
{
    std::vector<double> a = {1.0, 1.0};
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0}; // size 1, below minimum of 2
    std::vector<double> d = {1.0, 1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
}

// d.size() < 3: RHS vector has only two elements (line 113-114)
TEST(TRIDISOLVE_TEST, SIZE_ERROR_RHS_TOO_SMALL)
{
    std::vector<double> a = {1.0, 1.0};
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0}; // size 2, below minimum of 3
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
}

// a.size() != b.size()-1: subdiagonal is too long relative to diagonal (line 115-120)
TEST(TRIDISOLVE_TEST, SIZE_MISMATCH_SUBDIAG_VS_DIAG)
{
    std::vector<double> a = {1.0, 1.0, 1.0}; // size 3, but b.size()-1 == 2
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
}

// c.size() != b.size()-1: superdiagonal is too long relative to diagonal (line 121-126)
TEST(TRIDISOLVE_TEST, SIZE_MISMATCH_SUPERDIAG_VS_DIAG)
{
    std::vector<double> a = {1.0, 1.0};
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0, 1.0}; // size 3, but b.size()-1 == 2
    std::vector<double> d = {1.0, 1.0, 1.0};
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
}

// b.size() != d.size(): RHS has more elements than the diagonal (line 127-129)
TEST(TRIDISOLVE_TEST, SIZE_MISMATCH_DIAG_VS_RHS)
{
    std::vector<double> a = {1.0, 1.0};
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0};
    std::vector<double> d = {1.0, 1.0, 1.0, 1.0}; // size 4, but b.size() == 3
    std::vector<double> x;
    EXPECT_THROW(spida::tridisolve(a, b, c, d, x), detail::Exception);
}

// --- checkXInterp(vector) range-check tests (P1 gap coverage) ---

// LinearInterp::eval(vector): xi[0] is below x[0] (interp.cpp:14-18)
TEST(LINEAR_INTERP_TEST, ERROR_VECTOR_XINTERP_BELOW_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    std::vector<double> xi = {-0.1, 1.0}; // xi[0] < x[0]
    EXPECT_THROW((void) interp.eval(xi), detail::Exception);
}

// LinearInterp::eval(vector): xi.back() is above x.back() (interp.cpp:19-23)
TEST(LINEAR_INTERP_TEST, ERROR_VECTOR_XINTERP_ABOVE_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::LinearInterp interp(x, y);
    std::vector<double> xi = {0.5, 2.1}; // xi.back() > x.back()
    EXPECT_THROW((void) interp.eval(xi), detail::Exception);
}

// SplineInterp::eval(vector): xi[0] is below x[0] (interp.cpp:14-18)
TEST(SPLINE_INTERP_TEST, ERROR_VECTOR_XINTERP_BELOW_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::SplineInterp interp(x, y);
    std::vector<double> xi = {-0.1, 1.0}; // xi[0] < x[0]
    EXPECT_THROW((void) interp.eval(xi), detail::Exception);
}

// SplineInterp::eval(vector): xi.back() is above x.back() (interp.cpp:19-23)
TEST(SPLINE_INTERP_TEST, ERROR_VECTOR_XINTERP_ABOVE_RANGE)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    spida::SplineInterp interp(x, y);
    std::vector<double> xi = {0.5, 2.1}; // xi.back() > x.back()
    EXPECT_THROW((void) interp.eval(xi), detail::Exception);
}

// ============================================================
//  P2: SplineInterp constructor size<3 guard (interp.cpp:151-152)
// ============================================================

TEST(SPLINE_INTERP_TEST, ERROR_ONLY_TWO_POINTS_THROWS)
{
    // Two points passes checkData() but fails the spline minimum-size guard
    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {0.0, 1.0};
    EXPECT_THROW(spida::SplineInterp interp(x, y), detail::Exception);
}

TEST(SPLINE_INTERP_TEST, ERROR_ONE_POINT_THROWS)
{
    // One point fails checkData() directly
    std::vector<double> x = {0.0};
    std::vector<double> y = {0.0};
    EXPECT_THROW(spida::SplineInterp interp(x, y), detail::Exception);
}

TEST(SPLINE_INTERP_TEST, ERROR_SIZE_MISMATCH_THROWS)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0};
    EXPECT_THROW(spida::SplineInterp interp(x, y), detail::Exception);
}

TEST(SPLINE_INTERP_TEST, THREE_POINTS_MINIMUM_ACCEPTED)
{
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 4.0};
    EXPECT_NO_THROW(spida::SplineInterp interp(x, y));
    spida::SplineInterp interp(x, y);
    EXPECT_NEAR(interp.eval(1.0), 1.0, 1e-12);
}
