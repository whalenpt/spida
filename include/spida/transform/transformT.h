
#ifndef SPIDA_TRANSFORMST_H_
#define SPIDA_TRANSFORMST_H_ 

#include <vector>
#include "spida/constants.h"

namespace spida{

// interface class
class TransformsT 
{
public:
    TransformsT() {};
    virtual void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out) = 0;
    virtual void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out) = 0;
    virtual ~TransformsT() {};
};

}


#endif


