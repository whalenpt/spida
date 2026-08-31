/*------------------------------------------------------------------------------
 *
 *    Description: Worker-consumption example. This is the template a future
 *    spida-console job-service worker (a small executable linking
 *    SPIDA::spida directly, per the SPIDA Console architecture proposal's
 *    §1/§11) would follow: build a SimulationConfig (in a real worker,
 *    parsed from a request body — see spida/config/simulationconfig.h),
 *    run it via spida::config::SimulationRun, wire cooperative cancellation
 *    to SIGINT instead of killing the process, and print structured
 *    progress instead of relying on spdlog console output. Gated behind
 *    the SPIDA_WORKER_EXAMPLE CMake option (off by default).
 *
------------------------------------------------------------------------------*/

#include <spida/config/simulationbuilder.h>
#include <spida/config/simulationconfig.h>
#include <spida/helper/constants.h>

#include <csignal>
#include <iostream>

namespace {
spida::BasePropagator* g_propagator = nullptr;

extern "C" void handleSigint(int)
{
    if (g_propagator != nullptr)
        g_propagator->requestCancel();
}
} // namespace

int main()
{
    // A real worker parses this from the request body instead of
    // hardcoding it — see spida::config::SimulationConfig for the schema
    // (JSON-serializable via nlohmann::json, matching domain.ts's
    // SimulationConfig in the SPIDA Console proposal).
    spida::config::SimulationConfig cfg;
    cfg.name = "worker-example";
    cfg.grid.n = 8192;
    cfg.grid.a = -spida::PI;
    cfg.grid.b = spida::PI;
    cfg.solver.t0 = 0.0;
    cfg.solver.tf = 1.0;
    cfg.solver.hInit = 0.5;
    cfg.reporting.stepsPerOutput1D = 5;

    spida::config::SimulationRun run(cfg);
    g_propagator = &run.propagator();
    std::signal(SIGINT, handleSigint);

    // Structured progress instead of spdlog console lines — what a worker's
    // WS/event-log bridge would consume (see ProgressSnapshot).
    run.propagator().setProgressObserver([](const spida::ProgressSnapshot& s) {
        std::cout << "{\"t\":" << s.t;
        if (s.tf.has_value())
            std::cout << ",\"tf\":" << *s.tf;
        std::cout << ",\"stepsTaken\":" << s.stepsTaken << "}\n";
    });

    const bool step_ok = run.run();
    if (!step_ok) {
        std::cerr << "worker-example: a solver step failed (runtime_exception-equivalent)\n";
        return 1;
    }

    switch (run.propagator().stopReason()) {
    case spida::StopReason::CancelRequested:
        std::cout << "worker-example: cancelled\n";
        return 0;
    case spida::StopReason::Diverged:
        std::cerr << "worker-example: solution diverged\n";
        return 1;
    case spida::StopReason::MaxReportsReached:
    case spida::StopReason::None:
        std::cout << "worker-example: completed\n";
        return 0;
    }
    return 0;
}
