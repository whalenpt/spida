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

## Addendum: Phase C (nls_r / bessel_root_r — first complex-valued model, first non-uniform grid)

A fourth model, `nls_r` (radial cubic NLS), wired end to end — chosen over `kdv_cv` specifically
because `nls_r` already has a home in the existing wire contract (`GridKind::bessel_root_r`,
already in the enum) while `kdv_cv` would need an entirely new `GridKind` (a complex-valued
uniform grid) that exists in neither the C++ enum nor the proposal's own `domain.ts`. Wiring
`nls_r` exercises two things at once: the first complex-valued physical field, and the first
non-uniform grid — stress-testing whether `simulationbuilder.h`'s factory actually generalizes,
per this ADR's own original "Negative / follow-up" note.

1. **`spida::models::NlsR`/`NlsRPropagator`** (`include/spida/models/nls.h`), promoted from
   `demos/NLSR.cpp` following the exact pattern ADR-0003 established for
   `burgers`/`kdv`/`ks` — a working, self-contained demo already existed for this equation.
   `modelParams.gamma` (Kerr coefficient, default `2.0`) and `modelParams.amplitude` (initial
   Gaussian peak, default `2.0`) are exposed; the demo's fixed Gaussian width stays fixed.
   Reports on the grid's own (non-uniform) `r`/`kr` coordinates directly — no mirroring or
   synthesized axis, so the proposal's "gridCoords never assumed uniform" principle is
   satisfied without any extra bookkeeping.

2. **`SimulationRun`'s grid changed from a member to a per-case local variable.** The previous
   `spida::UniformGridRVX m_grid` member was constructed unconditionally, before `cfg.model`
   was even switched on — workable only because every wired model used the same grid type.
   Since `Burgers`/`Kdv`/`Ks`/`NlsR` all copy the grid into their own internal storage (see
   e.g. `kdv.h`'s `m_grid(grid)`), nothing needs `SimulationRun`'s own grid to outlive
   construction — so each `switch` case now declares its own concrete grid type
   (`UniformGridRVX` or `BesselRootGridR`) as a local, with no grid-side `std::variant` needed.
   `m_model` (which genuinely must outlive `m_propagator`/`m_solver`) stays a member variant,
   now including `spida::models::NlsR`.

3. **`validate()` (`include/spida/config/validation.h`) now checks model/grid *pairing*, not
   just whether each is independently wired** — directly implementing the proposal's own
   `config_validation` example ("incompatible model/grid pairing"), previously unreachable
   since only one `GridKind` was wired at all. `GridConfig.rMax` added (defaulted `5.0`,
   `BesselRootGridR`-shaped) and validated (`> 0`) the same way `grid.a < grid.b` is for
   uniform grids.

4. **`SimulationRun::propagator()`'s return type widened** from `BasePropagator&` to
   `PropagatorCV&` — source-compatible for every existing caller (all use only
   `BasePropagator`-declared members), and is what let a real numerical-verification test read
   the raw complex spectral field directly instead of needing a downcast.

5. **Numerical verification, done the way ADR-0003 set the bar — measured, not assumed.**
   `nls_r` has no closed-form solution for a general Gaussian, so it's checked against a
   conservation law: cubic NLS with a purely dispersive linear operator conserves the
   spectral-space L2 norm `sum_k |A_k|^2` exactly for the *continuous* equation. The first
   version of this check used a 1e-6 relative tolerance and **failed for real** — ~3e-4
   relative deviation at `tf=0.05`. Before loosening the tolerance blindly, this was
   investigated empirically: the deviation was unchanged between `grid.n=64` and `grid.n=256`
   (ruling out a spatial-resolution cause), and shrank roughly quadratically as `tf` shrank
   (~6e-6 at `tf=0.01`, ~1.4e-6 at `tf=0.005`) — confirming a genuine, shrinking
   time-integration truncation effect from the finite-order adaptive ETD scheme, not a fixed
   transform-normalization bug. Settled on `tf=0.01`, `epsRel=1e-10`, `1e-4` relative
   tolerance (~16x margin over the measured ~6e-6). See `test/config_tests.cpp`'s
   `NLS_R_POWER_IS_CONSERVED` for the full investigation trail in comments.

Verified via the same mandated Release build+test cycle as the addenda above, 11/11 binaries
passing, all new/updated tests passing individually.

## Consequences (Phase C)

- **Positive:** the `PropagatorCV`/variant factory pattern genuinely generalizes to a
  complex-valued field on a non-uniform grid without structural rework beyond the
  member-to-local grid change above — no grid-side variant, no changes to `BasePropagator`,
  `ReportHandler`, or the binary-frame/manifest work from Phase A.
- **Negative / follow-up:** `kdv_cv`/`nls_rt` remain unimplemented. `nls_rt` (2D, radial +
  time/frequency) would additionally need `Report2D`/`ReportComplex2D` wired through the
  worker's manifest — `docs/api/binary-frame-spec.md` already flags 2D framing as "not yet
  specified." `kdv_cv` still has no home in `GridKind` at all; wiring it means adding a new
  grid kind to both the C++ enum and the proposal's own `domain.ts`, a larger wire-contract
  change than this addendum took on.
