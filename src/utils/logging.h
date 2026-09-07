#pragma once

// spdlog is unavailable on the native, Conan-less Windows CI job (see
// .github/workflows/cmake.yml's "Windows (native — submodules, no Conan)"
// job — it configures straight cmake/Ninja with no `conan install` step, so
// find_package(spdlog CONFIG QUIET) in the top-level CMakeLists.txt never
// finds it and `if(TARGET spdlog::spdlog)` never links it into `spida`).
// Every other build (local dev per CLAUDE.md's mandatory sequence,
// Linux/macOS CI) goes through Conan, where spdlog is a hard requirement
// (conanfile.py) and this branch is always taken. __has_include lets this
// header compile either way instead of hard-failing on the Windows job,
// falling back to a no-op logger with the same ->info()/->debug()/->warn()/
// ->error() call surface so call sites never need to care which build
// they're in.
#if __has_include(<spdlog/spdlog.h>)
#define SPIDA_HAS_SPDLOG 1
#include <memory>
#include <spdlog/spdlog.h>
#else
#define SPIDA_HAS_SPDLOG 0
#endif

namespace detail {

#if SPIDA_HAS_SPDLOG

/// Lazily-created internal logger for spida's library-side diagnostics.
/// Writes to stderr (never stdout — some consumers, e.g. spida-worker,
/// treat stdout as a data channel). Not part of the public API; not
/// installed.
std::shared_ptr<spdlog::logger>& spidaLogger();

#else

/// No-op stand-in used when spdlog isn't available in this build (see the
/// __has_include guard above). Same ->info()/->debug()/->warn()/->error()
/// call surface as spdlog::logger.
struct NullLogger {
    template <typename... Args>
    void info(Args&&...) const noexcept
    {
    }
    template <typename... Args>
    void debug(Args&&...) const noexcept
    {
    }
    template <typename... Args>
    void warn(Args&&...) const noexcept
    {
    }
    template <typename... Args>
    void error(Args&&...) const noexcept
    {
    }
};

inline NullLogger* spidaLogger()
{
    static NullLogger logger;
    return &logger;
}

#endif

} // namespace detail
