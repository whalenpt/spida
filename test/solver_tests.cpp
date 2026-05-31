
#include <cmath>
#include <complex>
#include <filesystem>
#include <functional>
#include <vector>

#include <gtest/gtest.h>
#include <pwutils/pwexcept.h>
#include <spida/helper/constants.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/rkstiff/ETDCS.h>
#include <spida/rkstiff/IFAS.h>

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
    NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) { out[0] = in[0] * in[0]; };
}

// Exact solution for Bernoulli: a=-1, b=1, u0=0.5
static double bernoulliExact(double t)
{
    return 1.0 / (std::exp(t) + 1.0);
}

// Pure linear ODE: du/dt = lambda*u, NL=0  (for constant-step solvers)
static void makeLinearSystem(dcmplx lambda, LinOp& L, NLfunc& NL)
{
    L = {lambda};
    NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) { out[0] = dcmplx(0.0); };
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
static double bernoulli2ndExact(double t)
{
    return 1.0 / (1.5 * std::exp(2.0 * t) + 0.5);
}

// --- ETD34 adaptive-step solver tests ---

TEST(ETD34_TEST, BERNOULLI_ODE)
{
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
    makeBernoulli2D(L, NL);

    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-7);

    std::vector<dcmplx> u = {0.5, 0.5};
    double tf = 0.5;
    bool ok = solver.evolve(u, 0.0, tf, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(tf), 1e-5);
    EXPECT_NEAR(u[1].real(), bernoulli2ndExact(tf), 1e-5);
}

TEST(ETD34_TEST, CURRENT_TIME_UPDATED)
{
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
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
    // ETD35 (5th-order) must be strictly more accurate than ETD34 (4th-order) at the same epsRel
    EXPECT_LT(err34, 1e-3);
    EXPECT_LT(err35, 1e-3);
    EXPECT_LT(err35, err34);
}

// --- ETD4 constant-step solver tests ---

TEST(ETD4_TEST, PURE_LINEAR_DECAY)
{
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD4 solver(L, NL);

    std::vector<dcmplx> u = {1.0};
    solver.step(u, 0.1);
    EXPECT_NEAR(u[0].real(), std::exp(-0.1), 1e-8);
}

TEST(ETD4_TEST, PURE_LINEAR_OSCILLATION)
{
    // u' = i*omega*u -> |u(t)| = 1 for all t
    LinOp L;
    NLfunc NL;
    makeLinearSystem(spida::ii * spida::PI, L, NL);
    spida::ETD4 solver(L, NL);

    std::vector<dcmplx> u = {1.0};
    bool ok = solver.evolve(u, 0.0, 2.0, 0.01);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(std::abs(u[0]), 1.0, 1e-8);
}

TEST(ETD4_TEST, BERNOULLI_ODE)
{
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
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
    // NL = 0.01*u^2 is a small perturbation; modulus deviation is O(0.01) over t=[0,1]
    EXPECT_NEAR(std::abs(u[0]), 1.0, 1e-2);
}

// --- IF45DP adaptive-step Dormand-Prince integrating-factor solver tests ---

TEST(IF45DP_TEST, BERNOULLI_ODE)
{
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD34 solver(L, NL);

    EXPECT_NEAR(solver.modeCutoff(), 0.01, 1e-15);
    EXPECT_NEAR(solver.contourRadius(), 1.0, 1e-15);
    EXPECT_EQ(solver.contourPoints(), 32u);
}

TEST(SOLVER_ATTR_TEST, ETD34_CONTOUR_SETTERS)
{
    LinOp L;
    NLfunc NL;
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
    LinOp L;
    NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD4 solver(L, NL);

    EXPECT_NEAR(solver.modeCutoff(), 0.01, 1e-15);
    EXPECT_NEAR(solver.contourRadius(), 1.0, 1e-15);
    EXPECT_EQ(solver.contourPoints(), 32);
}

TEST(SOLVER_ATTR_TEST, SOLVER_SIZE)
{
    LinOp L = {-1.0, -2.0, -3.0};
    NLfunc NL = [](const std::vector<dcmplx>&, std::vector<dcmplx>& out) { out.assign(3, 0.0); };
    spida::ETD34 solver(L, NL);
    EXPECT_EQ(solver.size(), 3u);
}

TEST(SOLVER_ATTR_TEST, LOG_PROGRESS_DEFAULT)
{
    LinOp L;
    NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_FALSE(solver.logProgress());
    solver.setLogProgress(true);
    EXPECT_TRUE(solver.logProgress());
}

// ============================================================
//  Control threshold / error validation tests
// ============================================================

TEST(SOLVER_CONTROL_TEST, SET_INCREMENT_THRESHOLD_TOO_HIGH_THROWS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_THROW(solver.setIncrementThreshold(5.0), pw::Exception); // > MAX_S=4
}

TEST(SOLVER_CONTROL_TEST, SET_INCREMENT_THRESHOLD_TOO_LOW_THROWS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_THROW(solver.setIncrementThreshold(0.5), pw::Exception); // < 1.0
}

TEST(SOLVER_CONTROL_TEST, SET_INCREMENT_THRESHOLD_VALID_NO_THROW)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_NO_THROW(solver.setIncrementThreshold(2.0));
}

TEST(SOLVER_CONTROL_TEST, SET_DECREMENT_THRESHOLD_TOO_LOW_THROWS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_THROW(solver.setDecrementThreshold(0.1), pw::Exception); // < MIN_S=0.25
}

TEST(SOLVER_CONTROL_TEST, SET_DECREMENT_THRESHOLD_TOO_HIGH_THROWS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_THROW(solver.setDecrementThreshold(1.0), pw::Exception); // >= 1.0
}

TEST(SOLVER_CONTROL_TEST, SET_DECREMENT_THRESHOLD_VALID_NO_THROW)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_NO_THROW(solver.setDecrementThreshold(0.5));
}

TEST(SOLVER_CONTROL_TEST, SET_EPS_REL_NEGATIVE_THROWS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_THROW(solver.setEpsRel(-1e-6), pw::Exception);
}

TEST(SOLVER_CONTROL_TEST, SET_EPS_REL_ZERO_NO_THROW)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_NO_THROW(solver.setEpsRel(0.0));
}

// ============================================================
//  Solver accessor tests
// ============================================================

TEST(SOLVER_ACCESSOR_TEST, L_ACCESSOR_SIZE)
{
    LinOp L = {dcmplx(-1.0), dcmplx(-2.0)};
    NLfunc NL = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = in[0] * in[0];
        out[1] = in[1] * in[1];
    };
    spida::ETD34 solver(L, NL);
    // Access L() through base-class ref to avoid shadowing by ETD34's private member L
    spida::SolverCV& base = solver;
    EXPECT_EQ(base.L().size(), 2u);
    EXPECT_EQ(base.L()[0], dcmplx(-1.0));
}

TEST(SOLVER_ACCESSOR_TEST, DT_LAST_AFTER_EVOLVE)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-6);
    std::vector<dcmplx> u = {0.5};
    solver.evolve(u, 0.0, 1.0, 0.1);
    EXPECT_GT(solver.dtLast(), 0.0);
}

TEST(SOLVER_ACCESSOR_TEST, ACCEPT_AND_GET_Y_AFTER_STEP)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-6);
    std::vector<dcmplx> u = {0.5};
    double h = 0.1, h_next;
    solver.step(u, h, h_next);
    EXPECT_TRUE(solver.accept());
    EXPECT_EQ(solver.getY().size(), 1u);
    EXPECT_EQ(solver.getErr().size(), 1u);
}

TEST(SOLVER_ACCESSOR_TEST, NUM_THREADS_DEFAULT_AND_SET)
{
    LinOp L;
    NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD34 solver(L, NL);
    EXPECT_EQ(solver.numThreads(), 1u);
    solver.setNumThreads(2);
    EXPECT_EQ(solver.numThreads(), 2u);
}

TEST(SOLVER_ACCESSOR_TEST, FILE_REPORT_STATS_CREATES_FILE)
{
    namespace fs = std::filesystem;
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-6);
    std::vector<dcmplx> u = {0.5};
    solver.evolve(u, 0.0, 0.5, 0.1);

    fs::path dir = fs::temp_directory_path() / "spida_solver_stat_test";
    fs::create_directories(dir);
    solver.fileReportStats(dir);
    EXPECT_TRUE(fs::exists(dir / "solver_stat.dat"));
    fs::remove_all(dir);
}

// ============================================================
//  Evolve with PropagatorCV
// ============================================================

// Minimal propagator for solver integration tests
class SolverTestPropagator : public spida::PropagatorCV {
public:
    explicit SolverTestPropagator(const std::filesystem::path& path, unsigned sz = 1)
        : PropagatorCV(path), m_field(sz, dcmplx(0.5))
    {
    }

    void updateFields(double) override {}

    std::vector<dcmplx>& propagator() override
    {
        return m_field;
    }

private:
    std::vector<dcmplx> m_field;
};

TEST(EVOLVE_WITH_PROPAGATOR_TEST, ETD34_BERNOULLI_WITH_PROPAGATOR)
{
    namespace fs = std::filesystem;
    fs::path dir = fs::temp_directory_path() / "spida_etd34_prop_test";
    fs::create_directories(dir);

    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-7);

    SolverTestPropagator prop(dir);
    double tf = 1.0;
    bool ok = solver.evolve(prop, 0.0, tf, 0.1);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(prop.propagator()[0].real(), bernoulliExact(tf), 1e-5);
    fs::remove_all(dir);
}

TEST(EVOLVE_WITH_PROPAGATOR_TEST, ETD4_LINEAR_WITH_PROPAGATOR)
{
    namespace fs = std::filesystem;
    fs::path dir = fs::temp_directory_path() / "spida_etd4_prop_test";
    fs::create_directories(dir);

    LinOp L;
    NLfunc NL;
    makeLinearSystem(-1.0, L, NL);
    spida::ETD4 solver(L, NL);

    SolverTestPropagator prop(dir);
    prop.propagator()[0] = dcmplx(1.0); // u0 = 1, exact: u(t) = exp(-t)
    double tf = 1.0;
    bool ok = solver.evolve(prop, 0.0, tf, 0.01);

    EXPECT_TRUE(ok);
    EXPECT_NEAR(prop.propagator()[0].real(), std::exp(-1.0), 1e-5);
    fs::remove_all(dir);
}

// ============================================================
//  Evolve returns false (step size collapses)
// ============================================================

TEST(EVOLVE_FAILURE_TEST, ETD34_RETURNS_FALSE_WHEN_STEP_COLLAPSES)
{
    // NL so large that error is always huge → step halves until h < MIN_H
    LinOp L = {dcmplx(-1.0)};
    NLfunc NL_huge = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = dcmplx(1e50) * in[0];
    };
    spida::ETD34 solver(L, NL_huge);
    solver.setEpsRel(1e-10);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.1);
    EXPECT_FALSE(ok);
}

TEST(EVOLVE_FAILURE_TEST, IF34_RETURNS_FALSE_WHEN_STEP_COLLAPSES)
{
    LinOp L = {dcmplx(-1.0)};
    NLfunc NL_huge = [](const std::vector<dcmplx>& in, std::vector<dcmplx>& out) {
        out[0] = dcmplx(1e50) * in[0];
    };
    spida::IF34 solver(L, NL_huge);
    solver.setEpsRel(1e-10);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 1.0, 0.1);
    EXPECT_FALSE(ok);
}

// ============================================================
//  Solver with logProgress enabled (exercises logging branches)
// ============================================================

TEST(SOLVER_LOG_TEST, ETD34_EVOLVE_WITH_LOG_PROGRESS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD34 solver(L, NL);
    solver.setEpsRel(1e-6);
    solver.setLogProgress(true);
    solver.setLogFrequency(5);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 0.5, 0.1);
    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(0.5), 1e-4);
}

TEST(SOLVER_LOG_TEST, ETD4_EVOLVE_WITH_LOG_PROGRESS)
{
    LinOp L;
    NLfunc NL;
    makeBernoulliSystem(L, NL);
    spida::ETD4 solver(L, NL);
    solver.setLogProgress(true);
    solver.setLogFrequency(10);

    std::vector<dcmplx> u = {0.5};
    bool ok = solver.evolve(u, 0.0, 0.5, 0.01);
    EXPECT_TRUE(ok);
    EXPECT_NEAR(u[0].real(), bernoulliExact(0.5), 1e-4);
}
