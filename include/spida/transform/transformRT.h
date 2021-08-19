
#ifndef SPIDA_TRANSFORMSRT_H_
#define SPIDA_TRANSFORMSRT_H_ 

#include <vector>
#include "spida/constants.h"

namespace spida{

// interface class
class TransformsRT 
{
public:
    TransformsRT() {};
    virtual void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) = 0;
    virtual void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) = 0;
    virtual ~TransformsRT() {};
};



}


#endif


