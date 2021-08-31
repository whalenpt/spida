
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
        std::vector<dcmplx>& spectralField() {return m_spectral_field;}
        void initFields(const ShapeT& shape);
        PeriodicBLT& spida(); 
    private:
        ModelPeriodicBLT& m_md;
        std::vector<dcmplx> m_spectral_field;
        std::vector<double> m_real_field;
};

}


#endif






