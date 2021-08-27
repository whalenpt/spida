
#ifndef PROPAGATORT_H_
#define PROPAGATORT_H_

#include "spida/propagator/propagator.h"

namespace spida{

class ModelPeriodicBLT;
class PeriodicBLT;
class ShapeT;

class PropagatorBLT : public PropagatorCV
{
    public:
        PropagatorBLT(ModelPeriodicBLT& md); 
        virtual ~PropagatorBLT() {} 
        virtual void updateFields(double z) = 0;
        std::vector<double>& realField() {return m_real_field;}
        std::vector<dcmplx>& propagator() {return m_propagator;}
        void initFields(const ShapeT& shape);
        spida::PeriodicBLT& spida() {return m_spi;}
    private:
        spida::PeriodicBLT& m_spi;
        std::vector<dcmplx> m_propagator;
        std::vector<double> m_real_field;
};

}


#endif






