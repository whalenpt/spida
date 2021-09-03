
#ifndef PROPAGATORT_H_
#define PROPAGATORT_H_

#include "spida/PeriodicBLT.h"
#include "spida/propagator/propagator.h"
#include "spida/shape/shapeT.h"
#include "spida/model/modelT.h"

namespace spida{

class PropagatorBLT : public PropagatorCV
{
    public:
        PropagatorBLT(ModelPeriodicBLT& md); 
        virtual ~PropagatorBLT() {} 
        virtual void updateFields(double z) = 0;
        std::vector<double>& realField() {return m_real_field;}
        std::vector<dcmplx>& spectralField() {return m_spectral_field;}
        void initFields(const ShapeT& shape);
        PeriodicBLT& spida() {return m_md.spida();} 
        std::vector<dcmplx>& propagator() {return m_spectral_field;}
    private:
        ModelPeriodicBLT& m_md;
        std::vector<dcmplx> m_spectral_field;
        std::vector<double> m_real_field;
};

}


#endif






