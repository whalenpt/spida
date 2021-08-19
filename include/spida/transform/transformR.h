
#ifndef SPIDA_TRANSFORMR_H_
#define SPIDA_TRANSFORMR_H_ 

#include <vector>
#include "spida/grid/gridR.h"

namespace spida{

// interface class
class TransformR 
{
public:
    TransformR(const GridR& grid) {}
    virtual ~TransformR() {}
    TransformR()=delete;
    TransformR(const TransformR&)=delete;
    TransformR& operator=(const TransformR&)=delete;
    virtual void R_To_SR(const std::vector<double>& in,std::vector<double>& out) = 0;
    virtual void SR_To_R(const std::vector<double>& in,std::vector<double>& out) = 0;
};



}


#endif


