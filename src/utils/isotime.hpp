#pragma once

// Single shared ISO-8601 UTC timestamp formatter. Originally hand-rolled
// twice: once in worker/src/main.cpp (SimulationEvent "at", status.json
// timestamps) and, as of the "generatedAt" per-report-frame metadata field
// (see src/utils/report.hpp's buildMeta()), a second time would have been
// needed at the library level too. Pulled out here instead so both call
// sites share one gmtime_r/gmtime_s portability fork rather than risking the
// two drifting apart the way modelregistry.h's own header comment warns
// three independent hand-copies of the same fact eventually do.

#include <chrono>
#include <cstdio>
#include <ctime>
#include <string>

namespace detail {

/// Current wall-clock time, formatted "YYYY-MM-DDTHH:MM:SSZ" (UTC) — matches
/// docs/api/events.schema.json's "at" and openapi.yaml's date-time fields.
[[nodiscard]] inline std::string nowIso8601Utc()
{
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    // gmtime_r is POSIX-only -- MinGW/MSVC provide the thread-safe
    // equivalent as gmtime_s instead, with the arguments reversed (tm*
    // first, time_t* second) and an errno_t return instead of a tm*
    // return. See worker/src/main.cpp's git history for the Windows CI
    // build this distinction was first caught by.
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

} // namespace detail
