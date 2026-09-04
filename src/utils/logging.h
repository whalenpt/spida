#pragma once

#include <memory>

#include <spdlog/spdlog.h>

namespace detail {

/// Lazily-created internal logger for spida's library-side diagnostics.
/// Writes to stderr (never stdout — some consumers, e.g. spida-worker,
/// treat stdout as a data channel). Not part of the public API; not
/// installed.
std::shared_ptr<spdlog::logger>& spidaLogger();

} // namespace detail
