//constants.h
#pragma once

#include <complex>
#include <limits>
#include <map>

namespace spida{
    constexpr auto PI = 3.141592653589793238462643383279502884197;
	using dcmplx = std::complex<double>;
    constexpr dcmplx ii (0.0,1.0);
    constexpr auto NEAR_ZERO = 1.0e2*std::numeric_limits<double>::min();
	using stringMap = std::map<std::string,std::string>;
	using stringPair = std::pair<std::string,std::string>;
    enum class Dimension{D1,D2,D3};
}


