#include "spida/config/capabilities.h"
#include "spida/config/simulationbuilder.h"
#include "spida/config/simulationconfig.h"
#include "spida/config/validation.h"

#include <filesystem>
#include <vector>

#include <gtest/gtest.h>

namespace fs = std::filesystem;
using namespace spida::config;

// ============================================================
//  SimulationConfig <-> JSON round trip
// ============================================================

TEST(SIMULATION_CONFIG_TEST, DEFAULT_ROUND_TRIPS_THROUGH_JSON)
{
    SimulationConfig cfg;
    nlohmann::json j = cfg;
    SimulationConfig parsed = j.get<SimulationConfig>();

    EXPECT_EQ(parsed.name, cfg.name);
    EXPECT_EQ(parsed.model, cfg.model);
    EXPECT_EQ(parsed.grid.kind, cfg.grid.kind);
    EXPECT_EQ(parsed.grid.n, cfg.grid.n);
    EXPECT_DOUBLE_EQ(parsed.solver.epsRel, cfg.solver.epsRel);
    EXPECT_EQ(parsed.solver.kind, cfg.solver.kind);
    EXPECT_EQ(parsed.modelParams, cfg.modelParams);
    EXPECT_EQ(parsed.reporting.maxReports1D, cfg.reporting.maxReports1D);
    EXPECT_EQ(parsed.reporting.logFrequency, cfg.reporting.logFrequency);
}

TEST(SIMULATION_CONFIG_TEST, PARSES_PARTIAL_JSON_WITH_DEFAULTS)
{
    // Only override a couple of fields — everything else should fall back
    // to SimulationConfig's own defaults (NLOHMANN_DEFINE_TYPE_..._WITH_DEFAULT).
    nlohmann::json j = {{"name", "sweep-01"}, {"modelParams", {{"mu", 0.001}}}};
    SimulationConfig cfg = j.get<SimulationConfig>();

    EXPECT_EQ(cfg.name, "sweep-01");
    EXPECT_DOUBLE_EQ(cfg.modelParams.value("mu", -1.0), 0.001);
    EXPECT_EQ(cfg.model, ModelKind::burgers); // default
    EXPECT_EQ(cfg.grid.n, 256u);              // default
}

TEST(SIMULATION_CONFIG_TEST, MODEL_KIND_SERIALIZES_TO_EXPECTED_STRING)
{
    nlohmann::json j = ModelKind::burgers;
    EXPECT_EQ(j.get<std::string>(), "burgers");
}

// ============================================================
//  validate() (Phase B: structured, field-keyed config validation)
// ============================================================

TEST(VALIDATE_TEST, DEFAULT_CONFIG_IS_VALID)
{
    SimulationConfig cfg;
    EXPECT_TRUE(validate(cfg).empty());
}

TEST(VALIDATE_TEST, REJECTS_UNWIRED_MODEL_KIND_WITH_FIELD)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_r;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "model");
}

TEST(VALIDATE_TEST, REJECTS_UNWIRED_GRID_KIND_WITH_FIELD)
{
    SimulationConfig cfg;
    cfg.grid.kind = GridKind::cheb_x;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.kind");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_GRID_N)
{
    SimulationConfig cfg;
    cfg.grid.n = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.n");
}

TEST(VALIDATE_TEST, REJECTS_INVERTED_GRID_BOUNDS)
{
    SimulationConfig cfg;
    cfg.grid.a = 1.0;
    cfg.grid.b = -1.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.b");
}

TEST(VALIDATE_TEST, REJECTS_NEGATIVE_EPS_REL)
{
    // Regression check: this is exactly the case that used to fall through
    // to spida-worker's generic "runtime_exception" catch instead of
    // "config_validation" (Control::setEpsRel throws detail::Exception, not
    // std::invalid_argument) — see validation.h's header comment.
    SimulationConfig cfg;
    cfg.solver.epsRel = -1.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "solver.epsRel");
}

TEST(VALIDATE_TEST, REJECTS_NON_POSITIVE_H_INIT)
{
    SimulationConfig cfg;
    cfg.solver.hInit = 0.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "solver.hInit");
}

TEST(VALIDATE_TEST, REJECTS_TF_NOT_AFTER_T0)
{
    SimulationConfig cfg;
    cfg.solver.t0 = 1.0;
    cfg.solver.tf = 1.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "solver.tf");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_STEPS_PER_OUTPUT_1D)
{
    SimulationConfig cfg;
    cfg.reporting.stepsPerOutput1D = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.stepsPerOutput1D");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_MAX_REPORTS_1D)
{
    SimulationConfig cfg;
    cfg.reporting.maxReports1D = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.maxReports1D");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_LOG_FREQUENCY_ONLY_WHEN_LOG_PROGRESS_ENABLED)
{
    SimulationConfig cfg;
    cfg.reporting.logFrequency = 0;
    EXPECT_TRUE(validate(cfg).empty()); // logProgress defaults to false

    cfg.solver.logProgress = true;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.logFrequency");
}

TEST(VALIDATE_TEST, REPORTS_MULTIPLE_ERRORS_AT_ONCE)
{
    SimulationConfig cfg;
    cfg.grid.n = 0;
    cfg.solver.hInit = -1.0;
    EXPECT_EQ(validate(cfg).size(), 2u);
}

TEST(VALIDATE_TEST, VALIDATION_ERROR_ROUND_TRIPS_THROUGH_JSON)
{
    ValidationError err{"grid.n", "grid point count must be greater than zero"};
    nlohmann::json j = err;
    EXPECT_EQ(j.at("field"), "grid.n");
    ValidationError parsed = j.get<ValidationError>();
    EXPECT_EQ(parsed.field, err.field);
    EXPECT_EQ(parsed.message, err.message);
}

TEST(SIMULATION_RUN_CONSTRUCTOR_TEST, INVALID_CONFIG_STILL_THROWS_INVALID_ARGUMENT)
{
    // SimulationRun's constructor calls validate() itself (see
    // simulationbuilder.h) -- same std::invalid_argument contract as
    // before the split, now driven by validate() instead of a duplicate
    // set of checks.
    SimulationConfig cfg;
    cfg.solver.epsRel = -1.0;
    EXPECT_THROW(SimulationRun run(cfg), std::invalid_argument);
}

// ============================================================
//  Capability discovery (Phase B)
// ============================================================

TEST(CAPABILITIES_TEST, DESCRIBES_ALL_THREE_WIRED_MODELS)
{
    auto caps = spida::config::describeCapabilities();
    ASSERT_TRUE(caps.contains("models"));
    std::vector<std::string> names;
    for (auto const& m : caps.at("models"))
        names.push_back(m.at("model").get<std::string>());
    EXPECT_EQ(names, (std::vector<std::string>{"burgers", "kdv_rv", "ks"}));
}

TEST(CAPABILITIES_TEST, LISTS_ALL_FOUR_SOLVER_KINDS)
{
    auto caps = spida::config::describeCapabilities();
    ASSERT_TRUE(caps.contains("solvers"));
    EXPECT_EQ(caps.at("solvers").size(), 4u);
}

// ============================================================
//  SimulationRun (config -> grid/model/propagator/solver -> evolve)
// ============================================================

class SimulationRunTest : public ::testing::Test {
protected:
    void TearDown() override
    {
        // SimulationRun writes report files under "sim_<name>" in the CWD;
        // clean up so repeated test runs don't accumulate directories.
        std::error_code ec;
        fs::remove_all(fs::path("sim_" + m_cfg.name), ec);
    }

    SimulationConfig m_cfg;
};

TEST_F(SimulationRunTest, REJECTS_UNWIRED_MODEL_KIND)
{
    m_cfg.model = ModelKind::kdv_cv;
    EXPECT_THROW(SimulationRun run(m_cfg), std::invalid_argument);
}

TEST_F(SimulationRunTest, KDV_RV_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_kdv_rv";
    m_cfg.model = ModelKind::kdv_rv;
    m_cfg.modelParams = {{"solitonSpeed", 1.0}};
    m_cfg.grid.n = 128;
    m_cfg.grid.a = -20.0;
    m_cfg.grid.b = 20.0;
    m_cfg.solver.epsRel = 1e-10;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, KS_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_ks";
    m_cfg.model = ModelKind::ks;
    m_cfg.grid.n = 64;
    m_cfg.grid.a = 0.0;
    m_cfg.grid.b = 100.530964914873; // 32*pi, the reference domain
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.5;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, REJECTS_UNWIRED_GRID_KIND)
{
    m_cfg.grid.kind = GridKind::cheb_x;
    EXPECT_THROW(SimulationRun run(m_cfg), std::invalid_argument);
}

TEST_F(SimulationRunTest, BURGERS_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_to_completion";
    m_cfg.grid.n = 64;
    m_cfg.grid.a = -spida::PI;
    m_cfg.grid.b = spida::PI;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;
    m_cfg.reporting.maxReports1D = 500;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, PROGRESS_SNAPSHOT_TF_IS_POPULATED)
{
    // SimulationRun's constructor should call setFinalTime(cfg.solver.tf)
    // so every ProgressSnapshot a caller observes has tf set, not just t —
    // previously left unset (see simulationbuilder.h's header comment).
    m_cfg.name = "run_progress_tf";
    m_cfg.grid.n = 64;
    m_cfg.grid.a = -spida::PI;
    m_cfg.grid.b = spida::PI;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    std::vector<spida::ProgressSnapshot> seen;
    run.propagator().setProgressObserver(
        [&seen](const spida::ProgressSnapshot& s) { seen.push_back(s); });

    EXPECT_TRUE(run.run());
    ASSERT_FALSE(seen.empty());
    ASSERT_TRUE(seen.front().tf.has_value());
    EXPECT_DOUBLE_EQ(*seen.front().tf, m_cfg.solver.tf);
}

TEST_F(SimulationRunTest, CANCEL_REQUESTED_BEFORE_RUN_STOPS_EARLY)
{
    m_cfg.name = "run_cancelled";
    m_cfg.grid.n = 64;
    m_cfg.grid.a = -spida::PI;
    m_cfg.grid.b = spida::PI;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 1.0; // long enough that cancellation, not tf, ends it
    m_cfg.solver.hInit = 0.001;
    m_cfg.reporting.stepsPerOutput1D = 1;
    m_cfg.reporting.maxReports1D = 500;

    SimulationRun run(m_cfg);
    run.propagator().requestCancel();
    EXPECT_TRUE(run.run()); // solver-level evolve() still returns true (no step failed)
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::CancelRequested);
}
