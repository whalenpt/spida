#include "utils/logging.h"

#include <spdlog/sinks/stdout_color_sinks.h>

namespace detail {

std::shared_ptr<spdlog::logger>& spidaLogger()
{
    static std::shared_ptr<spdlog::logger> logger = [] {
        auto existing = spdlog::get("spida");
        return existing ? existing : spdlog::stderr_color_mt("spida");
    }();
    return logger;
}

} // namespace detail
