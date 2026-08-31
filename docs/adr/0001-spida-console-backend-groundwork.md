# ADR-0001: Backend groundwork for the SPIDA Console job-service proposal

## Status

Accepted — implemented and merged (PR #45, `feature/spida-console-backend-groundwork` → `main`, commit `c06f132`).

## Context

An architecture proposal ("SPIDA Console") describes a TypeScript frontend and job-service
backend for configuring, running, monitoring, and comparing SPIDA spectral-PDE simulations
from a browser. It is explicit that none of it exists today: SPIDA is a synchronous C++23
static library — build a grid/model/solver, call `evolve(...)` once, blocking, in-process,
writing report files to disk. There is no HTTP/RPC layer, no job/queue concept, no
cancellation primitive, and no structured status channel.

The proposal's own phased build order (§12) starts with "Phase 0 · Contract freeze" and
"Phase 1 · Prove the loop" in a **new, separate `spida-console` repo**, which would link
`SPIDA::spida` as an installed Conan package. Before any of that can be built, the library
itself needs a small number of extension seams — hooks a future worker process can use —
that do not exist in SPIDA today.

This ADR covers the groundwork done *inside the spida repo itself*, scoped deliberately to
stay within CLAUDE.md's constraints (static C++23 library, no service/HTTP layer, additive
changes only) while directly enabling the proposal's future worker.

## Decision

Add five additive, backward-compatible extension points, verified via the project's
mandated Release build+test cycle, and nothing else (no REST/WS/queue code, which belongs
in the separate `spida-console` repo per the proposal's own file structure):

1. **Cooperative cancellation** — `BasePropagator::requestCancel()` / `cancelRequested()`
   (atomic bool), checked at the existing `stepUpdate()` checkpoint that
   `SolverCV_AS`/`SolverCV_CS::evolve()` already gate their loop on. No solver-loop
   restructuring was needed — the checkpoint already existed, it just had nothing to check.

2. **Structured progress + divergence detection** — `ProgressSnapshot{t, tf, stepsTaken,
   currentStepSize}` pushed to an optional `setProgressObserver()` callback each step
   (`stepUpdate()` gained a `(t, dt)` overload so step size is available); a virtual
   `checkDiverged(t)` hook (default `false`) a subclass can override to flag a non-finite
   result. All three stop conditions (cancel / diverge / max-reports-reached) unify into one
   `StopReason` enum, queryable via `stopReason()`.

   **Not exceptions.** `SolverCV::evolve()` is `noexcept` end-to-end (confirmed by reading
   `solver.h`/`solver.cpp` before writing any code) — a thrown `DivergenceError`, the
   original idea, would call `std::terminate` rather than propagate. `StopReason` reuses the
   existing bool-return checkpoint contract instead.

3. **Reporting sink** — `report.hpp`'s report classes (`Report1D`, `ReportComplex1D`,
   `Report2D`, `ReportComplex2D`, `Track`, `TrackComplex`) now build a `nlohmann::json`
   (`buildJson()`) instead of writing directly to a stream; `ReportHandler::setSink(...)` /
   `setWriteReportFiles(false)` let a caller observe or replace the filesystem write with
   the same JSON object, exposed on `BasePropagator` as `setReportSink()` /
   `setWriteReportFiles()`. `nlohmann::json` was already a public, transitive dependency
   (`report.hpp` already serialized to JSON internally), so this added no new dependency.

4. **Config-driven construction** — `spida::config::SimulationConfig` (JSON-serializable,
   mirroring the proposal's `domain.ts` discriminated unions for `ModelKind`/`GridKind`/
   `SolverKind`) and `spida::config::SimulationRun`, which builds a grid/model/solver from
   it and exposes the same `evolve()` call every demo makes by hand today.

   **Scoped to a Burgers pilot only**, not all six models in the proposal's `ModelKind`
   enum — each model's math (`L`/`NL` operators) is genuinely hand-written C++, not data,
   and the proposal's own §12 Phase 1 recommends exactly this scope ("prove the loop... one
   model only — Burgers or KdV, closest to an existing demo"). `spida::models::Burgers` is a
   reusable, parameterized promotion of `demos/burgers.cpp`'s hand-coded classes; that demo
   was left untouched to avoid regression risk to a working, tested file. Extending to
   another model means adding its params struct to `simulationconfig.h`, a
   `spida::models::` type following the same pattern, and one `case` in the factory.

5. **`demos/worker_example.cpp`** (opt-in via new `SPIDA_WORKER_EXAMPLE` CMake option,
   default `OFF`) — exercises 1–4 together: config → `SimulationRun` → SIGINT-driven
   cooperative cancel → structured progress lines instead of spdlog console output. Placed
   under `demos/` rather than `usage/` specifically so it is covered by the same
   `-DSPIDA_TEST=ON -DSPIDA_DEMOS=ON` preset build CLAUDE.md mandates — `usage/` is a
   separate, non-preset standalone example project outside that build.

## Consequences

- **Positive:** a future `spida-worker` (in the separate `spida-console` repo) can parse a
  request into `SimulationConfig`, drive `SimulationRun`, cooperatively cancel it, and
  observe structured progress/report events — without SPIDA growing any transport
  dependency. 21 new GoogleTest cases (`propagator_tests.cpp`, new `config_tests.cpp`)
  cover the new seams; the existing test/demo suite was unaffected (verified via a clean
  `rm -rf build/` → `conan install --build=missing` → `cmake --preset conan-release` →
  build → `ctest` cycle, 11/11 binaries passing).
- **Negative / follow-up:** `StopReason` and `checkDiverged` only detect divergence a
  subclass explicitly checks for (the Burgers pilot checks `std::isfinite` on its physical
  field); there is no generic, library-level divergence detector. `SimulationConfig`/
  `SimulationRun` reject every `ModelKind` except `burgers` and every `GridKind` except
  `uniform_rvx` with `std::invalid_argument` — by design, but real work is deferred to
  whoever extends the pilot to KdV/KS/NLS.
- **Explicitly out of scope, on purpose:** no REST/WS/queue code was added anywhere; that
  work is [ADR-0002](0002-spida-console-phase2-live-feedback.md).

## Addendum: Phase B (structured validation, capability discovery)

Three follow-ups, done on the same branch as
[ADR-0003](0003-worker-relocation-and-cooperative-cancellation.md)'s addendum, closing gaps
identified by reading `SimulationRun`/`spida-worker` against this ADR's own stated goal —
"inline field errors on the config form; no run is created" — which plain exceptions can't
actually deliver:

1. **`spida::config::validate()`** (`include/spida/config/validation.h`), a pure function
   returning `std::vector<ValidationError>` (`{field, message}`, JSON-serializable). Mirrors
   every check `SimulationRun`'s constructor and the classes it builds (`UniformGridX`,
   `Control::setEpsRel`, `BasePropagator::validatePositive()`) already enforced via
   exceptions — but as data, keyed by field, instead of one free-text message. `SimulationRun`
   now calls it itself first (a private delegating constructor evaluates
   `requireValid(cfg)` before `m_grid` is even constructed), so it stays the single source of
   truth rather than a second, possibly-drifting opinion layered on top — still throws
   `std::invalid_argument` on failure, same contract as before.

   **A real, previously-silent bug this closes:** `Control::setEpsRel` throws
   `detail::Exception` (`src/utils/except.h`), not `std::invalid_argument` — `spida-worker`'s
   catch chain only mapped `std::invalid_argument`/`nlohmann::json::exception` to
   `failureReason: "config_validation"`, so a negative `epsRel` fell through to the generic
   `"runtime_exception"` catch instead, misclassifying a plain bad-input error as a crash.
   `validate()` rejects it before `Control::setEpsRel` is ever reached. Covered by
   `VALIDATE_TEST.REJECTS_NEGATIVE_EPS_REL`'s regression comment in `test/config_tests.cpp`.

   `spida-worker` now calls `validate()` immediately after parsing `config.json`, before
   writing `status: "running"` or constructing anything, and writes a structured
   `"validationErrors"` array into `status.json` alongside the existing `"failureDetail"`
   string (`docs/api/openapi.yaml`'s `Error`/`Simulation` schemas updated to match).

2. **`spida::config::describeCapabilities()`** (`include/spida/config/capabilities.h`),
   returning JSON: which `ModelKind`s are wired, each one's supported `GridKind`s and
   `modelParams` schema, and the full `SolverKind` list. Wired into `spida-worker --describe`
   (prints the JSON, exits 0, no config/output-dir required). Motivated directly by
   ADR-0003's own account of real drift: *"the wire shape had to be realigned... rather than
   the shape ADR-0001 had speculatively designed."* Lets a frontend/api-server introspect a
   worker binary instead of hand-syncing enums across the repo boundary — `docs/api/
   openapi.yaml` documents a (not-yet-implemented) `GET /capabilities` proposed to be backed
   by exactly this. **Not generated** — kept by hand in sync with `validate()` and
   `SimulationRun`'s own switch statements; extend it in the same commit that wires a new
   `ModelKind`/`GridKind`.

3. **`ReportingConfig`'s wire shape frozen now**, not extended later:
   `stepsPerOutput2D`/`maxReports2D`/`stepsPerOutputTrack` added (defaulted, currently unused
   — no wired model reports 2D/track data), matching the proposal's own `domain.ts` shape in
   full. Done specifically to avoid a second wire-shape break of the kind ADR-0003 already
   describes for `modelParams`/`logFrequency`'s placement.

Verified via the project's mandated Release build+test cycle (`conan install` → `cmake
--preset conan-release -DSPIDA_TEST=ON -DSPIDA_DEMOS=ON -DSPIDA_WORKER=ON` → build → `ctest`).
