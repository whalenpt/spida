# ADR-0003: Worker relocation into spida, and real cooperative cancellation

## Status

Accepted — implemented on `feature/spida-console-backend-phase2`.

## Context

Per [ADR-0002](0002-spida-console-phase2-live-feedback.md)'s open question
("where does the actual service code live?"), the decision landed on: the
worker (a C++ executable linking `SPIDA::spida` directly) lives in *this*
repo; the API server + job queue + metadata store (REST/WS, not C++) stays
in `spida-console` alongside the frontend. Rationale: the worker changes in
lockstep with spida's own API, needs no toolchain beyond what this repo
already has, and stays a plain CLI executable — no service/HTTP layer added
here, consistent with CLAUDE.md's scope.

A real worker already existed, authored directly in spida-console's
`services/worker` before this decision was made, and staged into this repo
at `temp/services/worker` for the move. It was not a prototype: three
models (Burgers, real-valued KdV, Kuramoto-Sivashinsky) were implemented
and numerically verified against closed-form solutions / linear stability
theory (see `worker/README.md`'s "Numerical verification" section for what
was checked), and it had been run for real against an installed spida
package, including catching and fixing two genuine packaging gaps in
spida's own public headers (`kiss_fftr.h` and an unused Boost include
leaking into every consumer's translation unit) along the way.

Critically, that worker's own code comments documented exactly the gaps
[ADR-0001](0001-spida-console-backend-groundwork.md) was written to close,
written *before* ADR-0001 existed:

- `src/main.cpp`'s header comment: *"the design doc imagined a custom
  ReportHandler hooking `BasePropagator::stepUpdate()` directly... That
  hook doesn't exist in spida 0.1.1's actual public API... the real fix
  belongs upstream, exposing `stepUpdate()`'s data the way the design doc
  assumed."*
- `jobs.ts`'s `cancelSimulation()` header comment: *"every
  `SolverCV::evolve()` overload is declared `noexcept`, so there is no safe
  way to unwind out of it early... What the design doc called
  'cooperative, checked at report checkpoints' isn't achievable without an
  upstream spida change."* Its workaround: a hard `SIGTERM`, a grace
  period, then `SIGKILL`.

ADR-0001's seams — `setProgressObserver()`, `requestCancel()`/
`cancelRequested()`, `StopReason`, `checkDiverged()` — are exactly that
upstream change, built independently and merged before this worker was
staged for relocation. This ADR is where the two meet.

## Decision

1. **Relocated** the worker to `worker/` (sibling to `test/`/`demos/`/
   `usage/`), as an in-tree CMake target (`add_executable(spida-worker
   src/main.cpp)`, `target_link_libraries(... PRIVATE SPIDA::spida)`),
   gated by a new `SPIDA_WORKER` option — dropping the Conan sub-project
   (`conanfile.py`, `CMakeUserPresets.json`, `conan-recipes/spida`) that
   existed solely to resolve spida as an external, installed package.
2. **Promoted** its three models into the library itself —
   `include/spida/models/{burgers,kdv,ks}.h` — following
   `spida::models::Burgers`'s existing pattern from ADR-0001 exactly,
   including a `checkDiverged()` override on each propagator. The worker's
   `main.cpp` no longer duplicates model definitions; it consumes the
   library's.
3. **Rewired construction through `spida::config::SimulationRun`**
   (ADR-0001's factory), generalized from Burgers-only to a
   `std::variant`-held model dispatching on `burgers`/`kdv_rv`/`ks`. As a
   side effect, all four `SolverKind`s are now honored via
   `config.json`'s `solver.kind` — the original worker always hardcoded
   `ETD35` regardless of what a request asked for.
4. **Replaced the progress-estimation workaround** (`stepsTaken =
   framesEmitted * stepsPerOutput1D`, an admitted undercount; `dtLast()`
   read directly off the solver) with the real
   `BasePropagator::setProgressObserver()` — exact `stepsTaken`, exact
   `currentStepSize`, no estimation. Throttled in `main.cpp` to fire at
   report cadence (`stepsPerOutput1D`), matching the original event rate —
   the observer itself fires every accepted solver step, which is much
   finer-grained and would otherwise multiply `events.ndjson`'s write
   volume several-fold for no benefit anything currently reads.
5. **Replaced the manifest-scanning `stopReason` heuristic** with the real
   `spida::StopReason`, string-mapped to stay wire-compatible with the
   original `"tf_reached"`/`"max_reports_reached"` values, adding
   `"cancel_requested"`/`"diverged"` for the two reasons that couldn't
   previously be distinguished. `api-server`'s `jobs.ts` already passes
   `status.json`'s `stopReason` through opaquely, so this needed no change
   there to be wire-compatible — see the consequence below for what
   *does*.
6. **`SIGTERM` now triggers real cooperative cancellation.** `api-server`
   already sends `SIGTERM` first (see `jobs.ts`'s `hardKill()`) — the
   worker now handles it by calling `requestCancel()` (an async-signal-safe
   atomic-bool store) instead of relying entirely on the grace-period
   `SIGKILL` escalation. Verified for real: a `tf=150` KS run, sent
   `SIGTERM` ~0.3s in, exited with code `0` and wrote its own
   `status.json` (`"status": "completed", "stopReason:
   "cancel_requested"`) — not killed.
7. **Fixed a real bug caught only by running the worker end to end**, not
   by the unit tests: `SimulationRun`'s constructor hardcoded report
   output to `"sim_" + cfg.name` (a path relative to CWD), ignoring
   whatever output directory the caller actually wanted. `configtest`'s
   existing tests never caught this because they never checked *where*
   `SimulationRun` wrote files, only that it ran. Fixed by adding an
   explicit `outDir` parameter (defaulting to the old behavior for
   demo/test convenience); `worker/main.cpp` now passes its real
   `<output-dir>` argument through.
8. **Fixed a second bug caught the same way**: `SimulationConfig`'s default
   `modelParams` member was declared `nlohmann::json modelParams{
   nlohmann::json::object()};` — brace-init, which nlohmann's
   initializer_list constructor treats as a one-element JSON *array*
   containing an empty object (`[{}]`), not the empty object itself. Any
   code path that never overrides `modelParams` (a default-constructed
   `SimulationConfig` run directly, not via JSON parsing) crashed calling
   `.value()` on it. Two of `configtest`'s own `SimulationRunTest` cases
   caught this immediately on the first real build+run. Fixed with direct
   construction (`= nlohmann::json::object()`) instead of brace-init.
9. **Realigned `SimulationConfig`'s wire shape** to match the real,
   already-working contract `apps/web`/this worker use, rather than the
   shape ADR-0001 had speculatively designed:
   - `SimulationConfig.burgersParams` (a typed, Burgers-only struct)
     replaced by `SimulationConfig.modelParams` (`nlohmann::json`, generic
     — `{"mu": ...}` / `{"solitonSpeed": ...}` / `{}` depending on
     `model`), matching `@spida-console/domain`'s actual field name and
     shape.
   - `SolverConfig.logFrequency` moved to `ReportingConfig.logFrequency`,
     matching where the real worker/frontend actually nest it.
   - `docs/api/openapi.yaml` updated to match.

## Consequence not addressed here (api-server, out of scope for this repo)

`services/api-server/src/jobs.ts`'s exit handler (`child.on("exit", ...)`)
special-cases `sim.status === "cancelling"`: it assumes a cancelling job is
always killed before it can write its own terminal `status.json`, and
resolves straight to `"cancelled"` from its own bookkeeping — see that
function's own header comment, quoted above, for why that assumption was
correct *at the time it was written*.

It no longer is. With real cooperative cancellation, a cancelled run now
usually reaches its next checkpoint and exits normally — code `0`, a real
`status.json` with `"status": "completed", "stopReason":
"cancel_requested"` — well before any `SIGKILL` grace period would fire.
`jobs.ts`'s current logic discards that real status.json unread and
substitutes its own synthesized `"cancelled"` instead.

This repo's change makes real cancellation possible; making `api-server`
actually use it — read the worker's own `status.json` first if it exists,
falling back to the kill-based `"cancelled"` only if the process was
actually killed before writing one — is spida-console's to make, not this
repo's. Flagged here so it isn't lost.

## Consequences (this repo)

- **Positive**: three real, numerically-verified models now live in the
  public library (`spida::models::{Burgers,Kdv,Ks}`), not duplicated in a
  downstream consumer. `spida::config::SimulationRun` covers all three end
  to end, and real cooperative cancellation works end to end for the first
  time — proven by an actual `SIGTERM` test, not just unit tests of the
  primitive in isolation.
- **Negative / follow-up**: `GridKind` still only wires `uniform_rvx` (all
  three models use it); `kdv_cv`/`nls_r`/`nls_rt` remain unimplemented,
  same as ADR-0001 left them. The `api-server` consequence above is real
  and unaddressed by this repo alone.

## Addendum: two gaps closed (Phase A)

Two gaps in the worker as relocated above — found by reading the code
against ADR-0001's own seams, not by anything failing — were closed on the
same branch:

1. **`ProgressSnapshot.tf` was always null.** Nothing between
   `spida::config::SimulationRun`'s constructor and the worker ever called
   `BasePropagator::setFinalTime()`, even though `cfg.solver.tf` is known
   the moment a config is parsed. Every progress event's `tf` field was
   silently `null` — `events.schema.json` and `EventLog::progress()` both
   correctly treat it as optional, so nothing broke loudly, but a
   `ProgressPanel`'s `t/tf` bar would have had nothing to divide by.
   Fixed: `SimulationRun`'s constructor now calls
   `m_propagator->setFinalTime(cfg.solver.tf)`. Covered by
   `SimulationRunTest.PROGRESS_SNAPSHOT_TF_IS_POPULATED`
   (`test/config_tests.cpp`).

2. **The manifest was built by re-reading files the worker had just
   written**, not by using the report sink ADR-0001 added for exactly this
   purpose. `buildManifest()` scanned `outDir` after the run, regex-matched
   `"<name>_<i>.json"` filenames to recover a series' name/index, and
   re-parsed each file's `"type"` field — the same file-tailing/re-parsing
   pattern ADR-0002 explicitly rejected for the *progress* channel ("added
   latency, races on partial writes, and re-parsing JSON the same process
   just serialized"). The same argument applies to the manifest and wasn't
   carried through when this worker was written. Fixed: `worker/src/main.cpp`
   now registers a `ManifestBuilder` via `BasePropagator::setReportSink()`,
   accumulating each series' `kind`/`valueType`/`frameCount` in-process, at
   the moment `report1D()`/`report2D()` call the sink — no directory scan, no
   regex, no re-parsing. Report files are still written to disk as before
   (`setWriteReportFiles` stays at its default `true`); only the manifest's
   *construction* changed.

   **Known limitation carried forward, now explicit in code**: the sink's
   JSON `"type"` field (`"xy"`/`"xy_complex"`/`"xyz"`/`"xyz_complex"`)
   cannot distinguish a `Track` series from a one-shot `Report1D` — both
   serialize `"xy"`. Moot today (none of `burgers`/`kdv_rv`/`ks` register a
   Track report), but the day one does, this needs either a separate sink
   per report kind or a richer `ReportHandler::Sink` signature that says
   which kind called it.

Verified via the project's mandated Release build+test cycle (`rm -rf
build/` not required — no CMakeLists.txt changes; `conan install` →
`cmake --preset conan-release -DSPIDA_TEST=ON -DSPIDA_DEMOS=ON
-DSPIDA_WORKER=ON` → build → `ctest`), 11/11 binaries passing.

## Addendum: `timeout` implemented (the error taxonomy's fourth reason)

The proposal's error taxonomy (§1) has five `FailureReason`s. Before this addendum, only
three were ever actually produced anywhere in this repo: `config_validation` (ADR-0001's Phase
B), `runtime_exception`, `divergence`. `timeout` ("optional operator wall-clock cap exceeded")
is now implemented too — `worker_crash` remains correctly out of scope (a crashed process
can't write its own `status.json`; that's the future api-server detecting a dead process).

**Deliberately not a `SimulationConfig` field.** The proposal frames `timeout` as an
*operator* cap, not a user-submitted parameter — so it's an optional 3rd `spida-worker` CLI
argument (`spida-worker config.json output-dir [timeout-seconds]`), not a new
`SimulationConfig`/wire-contract field. This is a worker/deployment concern parallel to how
`api-server` already supplies `SIGTERM` externally (see this ADR's own Decision §6) — the
submitted config never carries either.

**Reuses the existing cooperative-cancellation mechanism, adds nothing to the library.** The
worker's own `setProgressObserver()` callback (already registered for `events.ndjson`
forwarding) now also checks elapsed wall-clock time on *every* invocation (not throttled by
`stepsPerOutput1D` the way event forwarding is — a coarse report cadence shouldn't delay
noticing a timeout) and calls `BasePropagator::requestCancel()` the same way `handleSigterm()`
already does. No new `spida::StopReason` value, no library change at all — the worker tracks
its *own* `timedOut` flag locally to know *why* `requestCancel()` fired, since the library-level
`StopReason::CancelRequested` doesn't (and shouldn't) distinguish "asked via SIGTERM" from
"asked because the worker's own clock ran out" — that distinction only matters for choosing
which `FailureReason` to report, a worker-level concern.

**Precedence**: checked after the existing `!stepOk`/`Diverged` checks (a genuine solver
failure or divergence is a more specific diagnosis than "also ran out of time"), before the
generic "completed" path every other stop reason falls through to.

Verified via the mandated Release build+test cycle, 11/11 binaries passing (library untouched,
so no test changes were needed there — `worker/src/main.cpp` has no unit-test harness of its
own; correctness here rests on the clean compile plus this addendum's own review of the
callback/precedence logic). A real end-to-end CLI run (a long `ks` simulation against a 1s
budget) was attempted for direct verification, matching this repo's own "verified for real, not
just asserted" precedent, but was declined by the operator's permission settings in this
session (as it was twice earlier, for the Phase A and Phase C worker CLI smoke tests) — flagged
here rather than silently skipped.
