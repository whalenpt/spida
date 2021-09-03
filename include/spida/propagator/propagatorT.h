
#ifndef PROPAGATORT_H_
#define PROPAGATORT_H_

#include "spida/propagator/propagator.h"
#include "spida/model/modelT.h"
#include "spida/PeriodicBLT.h"
#include "spida/shape/shapeT.h"

namespace spida{

class PropagatorBLT : public PropagatorCV
{
    public:
        PropagatorBLT(ModelPeriodicBLT* md); 
        virtual ~PropagatorBLT() {} 
        virtual void updateFields(double z) = 0;
        virtual std::vector<dcmplx>& propagator() = 0;
};

}


#endif






