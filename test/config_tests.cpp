#include "spida/config/simulationbuilder.h"
#include "spida/config/simulationconfig.h"

#include <filesystem>

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
