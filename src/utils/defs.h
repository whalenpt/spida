#pragma once

#include <complex>
#include <map>

namespace detail {
using dcmplx = std::complex<double>;
inline constexpr dcmplx ii{0.0, 1.0};
using metadataMap = std::map<std::string, std::string, std::less<>>;
using metadataPair = std::pair<std::string, std::string>;

inline constexpr auto XLABEL = "x";
inline constexpr auto YLABEL = "y";
inline constexpr auto ZLABEL = "z";
} // namespace detail
