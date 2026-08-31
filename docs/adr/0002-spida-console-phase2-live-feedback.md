# ADR-0002: SPIDA Console Phase 2 — Live feedback (WS channel wired to report checkpoints)

## Status

Proposed — not implemented. Two open decisions block a concrete task breakdown (see
"Open questions" below); this ADR records the shape of the problem and a recommendation,
not a committed design.

## Context

The SPIDA Console proposal's §12 phased build order names Phase 2 as:

> **Phase 2 · Live feedback** — WS channel wired to the report checkpoints; RunTray,
> ProgressPanel, LogsTab.

This depends on Phase 0 (contract freeze: event schema, binary frame spec) and Phase 1
(prove the loop: a backend skeleton that validates → submits → runs a job **in a worker
thread** → reports completed/failed only) — per the proposal's own text, Phase 1 runs the
worker as an in-process thread inside the API server, not a separate OS process. Neither
Phase 0 nor Phase 1 exists yet: there is no `spida-console` repo, no `services/` directory,
and no WebSocket code anywhere in this workspace (confirmed by search before writing this
ADR).

[ADR-0001](0001-spida-console-backend-groundwork.md) added extension seams to the spida
C++ library itself — `requestCancel()`, `setProgressObserver()`, `checkDiverged()`,
`StopReason`, `setReportSink()` — but nothing that emits over a network transport. Phase 2
is the layer that would consume those seams and expose them live to a browser.

## Decision (proposed)

Split Phase 2 into three sub-parts:

**2a — Worker-side event emission.** Because Phase 1 runs the worker as an in-process
thread, wiring "report checkpoints" to a live channel is an in-process callback, not IPC:
the worker thread's `BasePropagator::setProgressObserver()` (from ADR-0001) pushes a
`SimulationEvent{kind: "progress"}` directly into the API server's per-simulation event
channel; existing spdlog-gated log points map to `{kind: "log"}`; state transitions
(`queued → running → completed/failed/cancelled`) map to `{kind: "status"}`.

**2b — API server: per-simulation event channel.** An in-memory ring buffer of recent
`SimulationEvent`s per running simulation, bounded (e.g. last N events), backing both:
- `WS /simulations/:id/stream` — accept a subscription, replay the buffer since a cursor
  on connect, then push live events.
- `GET /simulations/:id/events?since=<cursor>` — the same buffer, for reconnect gap-fill
  (per the proposal's §7 "reconnect without gaps").

**2c — Frontend: RunTray, ProgressPanel, LogsTab.** One multiplexed WS connection per the
proposal's §7 (`{subscribe: id}`/`{unsubscribe: id}` messages, demultiplexed by
`simulationId`), a Zustand store for the live/ephemeral state (§8), reconciled into the
TanStack Query cache for the rest of the run's metadata.

## Open questions (blocking a concrete task breakdown)

1. **Where does this get built?** No `spida-console`/`services/` location exists in this
   workspace. Needs a repo (new, or a path within an existing one) before 2a–2c can be
   scaffolded.
2. **Does Phase 0/1 get stood up first**, or is Phase 2 being planned in the abstract for
   later? 2a assumes a Phase 1 worker-thread skeleton that doesn't exist yet.

## Consequences (if built as proposed)

- **Positive:** no IPC/queue transport needed for the worker → API-server leg, since Phase
  1's worker is an in-process thread — 2a is a direct function call, not a message broker.
  Reuses ADR-0001's seams as-is; no further changes to the spida C++ library are implied by
  this ADR.
- **Negative / follow-up:** the moment a worker becomes a separate OS process (implied by
  the original proposal's own diagram in §1, and likely needed once cancellation escalates
  to a hard kill on timeout — see the proposal's error taxonomy), 2a's in-process callback
  stops working and needs a real IPC transport. That redesign is out of scope for this ADR
  and would be its own follow-up decision.
- **Alternative considered and rejected for 2a:** the worker thread tailing/polling the
  JSON files `ReportHandler` already writes, instead of using ADR-0001's observer callback.
  Rejected as an inferior default — added latency, races on partial writes, and re-parsing
  JSON the same process just serialized — but noted here since it requires zero further
  spida changes if that trade-off is ever preferred.
