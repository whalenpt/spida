# ADR-0001: Library extension seams for an external job-service worker

## Status

Accepted — implemented.

## Context

SPIDA is a synchronous, blocking, in-process C++23 library: build a grid/model/solver
by hand, call `evolve(...)` once, block until it returns, report files land on disk.
An external job-service worker (`worker/`, see ADR-0002) that runs simulations on
behalf of a browser frontend needs several things this library didn't originally
expose: a way to cancel a run in progress, a way to observe progress as it happens
(not just after it finishes), a way to intercept report data without a filesystem
round-trip, and a way to build a grid/model/solver from a JSON config instead of
hand-written code.

## Decision

Five additive, backward-compatible extension points:

1. **Cooperative cancellation** — `BasePropagator::requestCancel()` /
   `cancelRequested()` (an atomic bool), checked at the existing `stepUpdate()`
   checkpoint `SolverCV_AS`/`SolverCV_CS::evolve()` already gate their loop on. No
   solver-loop restructuring needed — the checkpoint already existed, it just had
   nothing to check. `requestCancel()` only stores to the atomic, which is
   async-signal-safe — safe to call from a `SIGTERM` handler (see ADR-0003).

2. **Structured progress + divergence detection** — `ProgressSnapshot{t, tf,
   stepsTaken, currentStepSize}` pushed to an optional `setProgressObserver()`
   callback each step; a virtual `checkDiverged(t)` hook (default `false`) a model
   subclass can override to flag a non-finite result. Cancel/diverge/max-reports-reached
   all unify into one `StopReason` enum, queryable via `stopReason()`. **Not
   exceptions**: `SolverCV::evolve()` is `noexcept` end-to-end, so a thrown error
   would call `std::terminate` rather than propagate — `StopReason` reuses the
   existing bool-return checkpoint contract instead.

3. **Reporting sink** — `report.hpp`'s report classes (`Report1D`, `ReportComplex1D`,
   `Report2D`, `ReportComplex2D`, `Track`, `TrackComplex`) build a `nlohmann::json`
   (`buildJson()`) instead of writing directly to a stream; `ReportHandler::setSink(...)`
   / `setWriteReportFiles(false)` let a caller observe or replace the filesystem write
   with the same JSON object, exposed on `BasePropagator` as `setReportSink()` /
   `setWriteReportFiles()`.

4. **Config-driven construction** — `spida::config::SimulationConfig`
   (JSON-serializable) and `spida::config::SimulationRun`, which builds a
   grid/model/solver from it and exposes the same `evolve()` call every demo makes by
   hand. See ADR-0002 for how this grew from an initial Burgers-only pilot to cover
   all six wired models via a shared per-model registry.

5. **`demos/worker_example.cpp`** (opt-in via `SPIDA_WORKER_EXAMPLE`, default `OFF`) —
   exercises 1–4 together: config → `SimulationRun` → SIGINT-driven cooperative cancel
   → structured progress lines. Lives under `demos/` specifically so it's covered by
   the same `-DSPIDA_TEST=ON -DSPIDA_DEMOS=ON` preset build CLAUDE.md mandates.

## Consequences

- **Positive**: `worker/` (ADR-0002) drives `SimulationRun`, cooperatively cancels it,
  and observes structured progress/report events — without SPIDA growing any
  transport dependency of its own. Covered by GoogleTest cases in
  `test/propagator_tests.cpp`/`test/config_tests.cpp`.
- **Negative / accepted limitation**: `StopReason`/`checkDiverged` only detect
  divergence a model subclass explicitly checks for (e.g. `std::isfinite` on its
  physical field) — there is no generic, library-level divergence detector.
- **Out of scope, on purpose**: no REST/WS/queue/HTTP-server code belongs in this
  repo — that's `spida-console`'s, a separate repo that doesn't exist in this
  workspace yet. See the new CLAUDE.md section for the full scope boundary.
