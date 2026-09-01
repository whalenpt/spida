# ADR-0003: Transport and live events

## Status

Accepted, implemented on the worker side. The api-server/frontend side is
`spida-console`'s to build (a separate repo, does not exist in this workspace) —
flagged explicitly below wherever this ADR's decisions imply work there.

## Context

An early version of this design assumed a job-service worker would run as an
**in-process thread** inside the api-server, so "wire progress to a live channel"
would be a direct function call — no transport needed at all. That assumption turned
out to be wrong: `spida-worker` (ADR-0002) is a real, separate OS process — one job,
exits. Cooperative cancellation already proves this works end to end: `api-server`
sends `SIGTERM`, the worker's handler calls `requestCancel()` (ADR-0001), and the
process exits normally (code `0`) with its own `status.json`
(`"status": "completed", "stopReason": "cancel_requested"`) — not killed. What was
never resolved is the *live update* side of that same subprocess reality: how does an
api-server actually get status/progress/results out of a worker process as it runs,
and how does it get them to a browser.

## Decision

**Worker → api-server: the worker's own stdout, not file polling or watching.**
`EventLog` (`worker/src/main.cpp`) builds one `SimulationEvent` JSON object per
status/progress/log/report change and fans it out to two independent `EventSink`s
from a single producer, so they can't drift apart:

- `FileEventSink` → `events.ndjson`, appended — durable, replayable after the process
  exits or after an api-server restart (`GET /simulations/:id/events?since=`).
- `StdoutEventSink` → the worker's own stdout, one JSON line per event, **explicitly
  flushed**. A parent that spawns `spida-worker` already holds a pipe to its stdout
  the instant it does — no polling `status.json`, no filesystem-watching
  `events.ndjson` for new lines, none of the partial-write races that come with
  tailing a file mid-append. The explicit flush matters: stdout is fully buffered
  (not line-buffered) once redirected to a pipe rather than a tty, so without it
  events would sit in the process's internal buffer and only arrive in a burst at
  exit — silently defeating the point.

**`manifest.json` is rewritten live, on every report frame, not just once at the
end.** Without this, `GET /simulations/:id/results` can't answer until a run is
completely finished, which undercuts the entire live-transport point — a caller
could see "73% done" but not what the field looks like at 73%. A new `"report"`
`SimulationEvent` kind (`{seriesName, frameIndex}`, 0-based, matching the real
on-disk `<seriesName>_<frameIndex>.json`) fires alongside each rewrite, so a live
subscriber can fetch a new frame the instant it exists instead of polling.

**api-server → browser: Server-Sent Events, not WebSocket.** The data flow here is
one-directional (server pushes `SimulationEvent`s; cancel is already, correctly, a
one-shot `POST /simulations/:id/cancel`, not a stream message) — SSE is a structural
match, not WebSocket's full duplex carried for no reason. Its native `Last-Event-ID`
reconnect maps directly onto the `since=<cursor>` replay shape already in
`events.schema.json`/`openapi.yaml`, without hand-rolling a WS reconnect protocol.
Recommended shape end to end: plain HTTP for commands/queries (`POST /simulations`,
`POST /simulations/:id/cancel`, `GET /simulations/:id`), this stdout NDJSON pipe for
worker→api-server, SSE for api-server→browser.

**Exit codes** (previously implicit, now documented — see `worker/README.md`): `0` on
any terminal `"completed"` state (`tf` reached, `max_reports_reached`, **or**
`cancel_requested` — a cooperatively-cancelled run still exits `0`); `1` on a
terminal `"failed"` state; `2` on a CLI usage error, with no `status.json` written at
all. A process that exits with neither is a real crash (`worker_crash` in the error
taxonomy) — the worker can never report that about itself; detecting it is the
caller's job.

## Consequences

- **Positive**: for a single api-server instance, the whole pipeline is
  strikingly thin — read a line from the worker's stdout pipe, forward it as one SSE
  frame, append it to `events.ndjson`. No message broker, no queue, no WS server.
- **Negative / follow-up, not addressed here**: this in-process "tail my own child's
  stdout" model doesn't survive an api-server scaled to multiple instances — a
  browser connected to instance B can't see a worker instance A spawned. Fix then
  (sticky-route by simulation id first; a real pub/sub layer only if that's
  insufficient) — solving a scaling problem that doesn't exist yet would be
  premature.
- **Negative / follow-up, api-server's to make**: a cooperatively-cancelled run now
  usually reaches its next checkpoint and exits normally with a real terminal
  `status.json`, well before any hard-kill grace period would fire. An exit handler
  written under the old in-process-thread assumption — one that assumes a
  "cancelling" job is always killed before it can write its own terminal status, and
  synthesizes `"cancelled"` from its own bookkeeping instead — needs to read that
  real `status.json` first, falling back to a synthesized `"cancelled"` only if the
  process was actually killed before writing one.
- **Known, accepted gap**: the worker's `SIGTERM` handler is installed only after
  `SimulationRun` finishes constructing. A `SIGTERM` arriving before that point (config
  parsing/validation, or construction itself — milliseconds, at the problem sizes
  this repo runs) gets the default disposition: immediate termination, no
  `status.json` written, rather than cooperative cancellation. Not closed here —
  doing so needs a pre-construction cancel flag checked at a couple of points, more
  complexity than the size of the actual risk currently justifies. Worth knowing if
  an api-server integration test sends cancel immediately after submit and
  occasionally sees a run vanish with no terminal status.
