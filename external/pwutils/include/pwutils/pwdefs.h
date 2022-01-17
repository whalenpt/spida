// pwdefs.h
#pragma once

#include <map>
#include <complex>

namespace pw{
    using dcmplx = std::complex<double>;
    const dcmplx ii (0.0,1.0);
    using metadataMap = std::map<std::string,std::string>;
    using metadataPair = std::pair<std::string,std::string>;
    enum class FileSignature{DAT,JSON,UNKNOWN};
    enum class DataSignature{XY,XCVY,XYZ,XYCVZ,UNKNOWN};
    enum class OperatorSignature{NONE,LOGX,LOGY,LOGXLOGY,LOGZ};

    constexpr auto XLABEL = "x";
    constexpr auto YLABEL = "y";
    constexpr auto ZLABEL = "z";
}


