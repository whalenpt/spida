
#ifndef PROPAGATORST_H_
#define PROPAGATORST_H_

#include "spida/model/modelT.h"
#include "spida/propagator/propagator.h"
#include "spida/shape/shapeT.h"
#include "spida/spidaT.h"

namespace spida{

class PropagatorsT : public PropagatorsDC
{
    public:
        PropagatorsT(ModelT& md) : m_spi(md.spida()), 
            m_U(m_spi.getNt(),0.0), m_V(m_spi.getNst(),0.0) {}
        virtual ~PropagatorsT() {} 
        virtual void updateFields(double z) = 0;
        std::vector<double>& real() {return m_U;}
        std::vector<dcmplx>& spec() {return m_V;}
        void initFields(const ShapeT& shape);
        spida::PeriodicT& spida() {return m_spi;}
    private:
        spida::PeriodicT& m_spi;
        std::vector<double> m_U;
        std::vector<dcmplx> m_V;
};

}


#endif






