#ifndef SPIDA_CONSTANTS_H
#define SPIDA_CONSTANTS_H

#include <complex>
#include <limits>
#include <map>

namespace spida{
    const double PI = 3.141592653589793238462643383279502884197;
	using dcmplx = std::complex<double>;
    const dcmplx ii (0.0,1.0);
    const double NEAR_ZERO = 1.0e2*std::numeric_limits<double>::min();
	using stringMap = std::map<std::string,std::string>;
	using stringPair = std::pair<std::string,std::string>;
}



#endif
