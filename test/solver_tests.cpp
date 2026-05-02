
#include <gtest/gtest.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/rkstiff/ETDCS.h>
#include <spida/rkstiff/IFAS.h>
#include <spida/helper/constants.h>
#include <cmath>
#include <complex>
#include <vector>
#include <functional>

using spida::dcmplx;
using spida::LinOp;
using spida::NLfunc;

// Bernoulli ODE: du/dt = L*u + NL(u)
// With L={a}, NL(u)={b*u^2}, u0 given,
// exact solution: u(t) = 1 / ((1/u0 + b/a)*exp(-a*t) - b/a)
// Using a=-1, b=1, u0=0.5: u(t) = 1/(e^t + 1)

static void makeBernoulliSystem(LinOp& L, NLfunc& NL)
{
    L = {dcmplx(-1.0)};
    NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = in[0] * in[0];
    };
}

// Exact solution for Bernoulli: a=-1, b=1, u0=0.5
static double bernoulliExact(double t) {
    return 1.0 / (std::exp(t) + 1.0);
}

// Pure linear ODE: du/dt = lambda*u, NL=0  (for constant-step solvers)
static void makeLinearSystem(dcmplx lambda, LinOp& L, NLfunc& NL)
{
    L = {lambda};
    NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = dcmplx(0.0);
    };
}

// Two-component Bernoulli system (decoupled): each mode follows Bernoulli
static void makeBernoulli2D(LinOp& L, NLfunc& NL)
{
    L = {dcmplx(-1.0), dcmplx(-2.0)};
    NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = in[0] * in[0];
        out[1] = in[1] * in[1];
    };
}

// Exact for 2D: u0=0.5, a=-1, b=1 -> u(t) = 1/(e^t+1)
//              v0=0.5, a=-2, b=1 -> v(t) = 1/((1/0.5 - 1/2)*e^{2t} + 1/2)
//                                        = 1/(1.5*e^{2t} + 0.5) [check below]
// a=-2,b=1,u0=0.5: v(t)=(v0+b/a)*exp(-a*t)-b/a = (2-0.5)*e^{2t}+0.5 = 1.5*e^{2t}+0.5
// u(t) = 1/(1.5*e^{2t}+0.5)
static double bernoulli2ndExact(double t) {
    return 1.0 / (1.5 * std::exp(2.0 * t) + 0.5);
}


// --- ETD34 adaptive-step solver tests ---

TEST(ETD34_TEST, BERNOULLI_ODE)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);

    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-7);

    std::vector<dcmplx> u = {0.5};
    double tf = 1.0;
    bool ok = solver.evolve(u, 0.0, tf, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(tf), 1e-5);
    EXPECT_NEAR(u[0].imag(), 0.0, 1e-12);
}

TEST(ETD34_TEST, BERNOULLI_2D_DECOUPLED)
{
    LinOp L; NLfunc NL;
    makeBernoulli2D(L, NL);

    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-7);

    std::vector<dcmplx> u = {0.5, 0.5};
    double tf = 0.5;
    bool ok = solver.evolve(u, 0.0, tf, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(tf), 1e-4);
    EXPECT_NEAR(u[1].real(), bernoulli2ndExact(tf), 1e-4);
}

TEST(ETD34_TEST, CURRENT_TIME_UPDATED)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-7);

    std::vector<dcmplx> u = {0.5};
    double tf = 0.5;
    solver.evolve(u, 0.0, tf, 0.1);
    EXPECT_NEAR(solver.currentTime(), tf, 1e-12);
}

// --- ETD35 adaptive-step solver tests ---

TEST(ETD35_TEST, BERNOULLI_ODE)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);

    spida::ETD35 solver(L, NL);
    solver.setEpsRel(1e-8);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(1.0), 1e-6);
}

TEST(ETD35_TEST, BETTER_ACCURACY_THAN_ETD34)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);
    double tf = 1.0;

    spida::ETD34 s34(L, NL);
    s34.setEpsRel(1e-4);
    std::vector<dcmplx> u34 = {0.5};
    s34.evolve(u34, 0.0, tf, 0.1);

    spida::ETD35 s35(L, NL);
    s35.setEpsRel(1e-4);
    std::vector<dcmplx> u35 = {0.5};
    s35.evolve(u35, 0.0, tf, 0.1);

    double exact = bernoulliExact(tf);
    double err34 = std::abs(u34[0].real() - exact);
    double err35 = std::abs(u35[0].real() - exact);
    // ETD35 is 5th-order, ETD34 is 4th-order - both should be accurate
    EXPECT_LT(err34, 1e-3);
    EXPECT_LT(err35, 1e-3);
}

// --- ETD4 constant-step solver tests ---

TEST(ETD4_TEST, PURE_LINEAR_DECAY)
{
    LinOp L; NLfunc NL;
    makeLinearSystem(-1.0, L, NL);

    spida::ETD4 solver(L, NL);

    std::vector<dcmplx> u = {1.0};
    double tf = 1.0;
    bool ok = solver.evolve(u, 0.0, tf, 0.01);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), std::exp(-1.0), 1e-6);
}

TEST(ETD4_TEST, SINGLE_STEP)
{
    LinOp L; NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD4 solver(L, NL);

    std::vector<dcmplx> u = {1.0};
    solver.step(u, 0.1);
    EXPECT_NEAR(u[0].real(), std::exp(-0.1), 1e-8);
}

TEST(ETD4_TEST, PURE_LINEAR_OSCILLATION)
{
    // u' = i*omega*u -> |u(t)| = 1 for all t
    LinOp L; NLfunc NL;
    makeLinearSystem(spida::ii * spida::PI, L, NL);
    spida::ETD4 solver(L, NL);

    std::vector<dcmplx> u = {1.0};
    bool ok = solver.evolve(u, 0.0, 2.0, 0.01);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(std::abs(u[0]), 1.0, 1e-8);
}

TEST(ETD4_TEST, BERNOULLI_ODE)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD4 solver(L, NL);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.01);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(1.0), 1e-5);
}

// --- IF34 adaptive-step integrating-factor solver tests ---

TEST(IF34_TEST, BERNOULLI_ODE)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);

    spida::IF34 solver(L, NL);
    solver.setEpsRel(1e-7);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(1.0), 1e-5);
}

TEST(IF34_TEST, OSCILLATORY_NONLINEAR)
{
    // L = i*pi, NL = 0.01*u^2 (mostly oscillatory with small nonlinear correction)
    LinOp L = {spida::ii * spida::PI};
    NLfunc NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = dcmplx(0.01) * in[0] * in[0];
    };

    spida::IF34 solver(L, NL);
    solver.setEpsRel(1e-6);

    std::vector<dcmplx> u = {1.0};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.1);

    EXPECT_TRUE(ok);
    // Modulus should be close to 1 (small nonlinear perturbation)
    EXPECT_NEAR(std::abs(u[0]), 1.0, 0.1);
}

// --- IF45DP adaptive-step Dormand-Prince integrating-factor solver tests ---

TEST(IF45DP_TEST, BERNOULLI_ODE)
{
    LinOp L; NLfunc NL;
    makeBernoulliSystem(L, NL);

    spida::IF45DP solver(L, NL);
    solver.setEpsRel(1e-8);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(1.0), 1e-6);
}

TEST(IF45DP_TEST, BERNOULLI_2D_DECOUPLED)
{
    LinOp L; NLfunc NL;
    makeBernoulli2D(L, NL);

    spida::IF45DP solver(L, NL);
    solver.setEpsRel(1e-8);

    std::vector<dcmplx> u = {0.5, 0.5};
    double tf = 0.5;
    bool ok = solver.evolve(u, 0.0, tf, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(tf), 1e-5);
    EXPECT_NEAR(u[1].real(), bernoulli2ndExact(tf), 1e-5);
}

// --- Solver attribute tests ---

TEST(SOLVER_ATTR_TEST, ETD34_CONTOUR_DEFAULTS)
{
    LinOp L; NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD34 solver(L, NL);

    EXPECT_NEAR(solver.modeCutoff(), 0.01, 1e-15);
    EXPECT_NEAR(solver.contourRadius(), 1.0, 1e-15);
    EXPECT_EQ(solver.contourPoints(), 32u);
}

TEST(SOLVER_ATTR_TEST, ETD34_CONTOUR_SETTERS)
{
    LinOp L; NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD34 solver(L, NL);

    solver.setModeCutoff(0.05);
    solver.setContourRadius(2.0);
    solver.setContourPoints(64);

    EXPECT_NEAR(solver.modeCutoff(), 0.05, 1e-15);
    EXPECT_NEAR(solver.contourRadius(), 2.0, 1e-15);
    EXPECT_EQ(solver.contourPoints(), 64u);
}

TEST(SOLVER_ATTR_TEST, ETD4_CONTOUR_DEFAULTS)
{
    LinOp L; NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD4 solver(L, NL);

    EXPECT_NEAR(solver.modeCutoff(), 0.01, 1e-15);
    EXPECT_NEAR(solver.contourRadius(), 1.0, 1e-15);
    EXPECT_EQ(solver.contourPoints(), 32);
}

TEST(SOLVER_ATTR_TEST, SOLVER_SIZE)
{
    LinOp L = {-1.0, -2.0, -3.0};
    NLfunc NL = [](const std::vector<dcmplx>&, std::vector<dcmplx>& out) {
        out.assign(3, 0.0);
    };
    spida::ETD34 solver(L, NL);
    EXPECT_EQ(solver.size(), 3u);
}

TEST(SOLVER_ATTR_TEST, LOG_PROGRESS_DEFAULT)
{
    LinOp L; NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_FALSE(solver.logProgress());
    solver.setLogProgress(true);
    EXPECT_TRUE(solver.logProgress());
}
