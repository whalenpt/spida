
#include <cmath>
#include <memory>
#include "spida/RCVT.h"
#include "spida/grid/uniformCVT.h"
#include "spida/grid/besselR.h"
#include "spida/transform/hankelfftRCVT.h"
#include <pwutils/pwexcept.h>

namespace spida{

SpidaRCVT::SpidaRCVT(const BesselRootGridR& gridR,\
        const UniformGridCVT& gridT,unsigned threads) :
    m_gridR(std::make_unique<BesselRootGridR>(gridR)),
    m_gridT(std::make_unique<UniformGridCVT>(gridT)),
    m_tr(std::make_unique<HankelFFTRCVT>(gridR,gridT,threads)) { }
  

const BesselRootGridR& SpidaRCVT::getGridR() const { return *m_gridR; }
const UniformGridCVT& SpidaRCVT::getGridT() const { return *m_gridT; }
const HankelFFTRCVT& SpidaRCVT::getTransformRT() const { return *m_tr; }
const std::vector<double>& SpidaRCVT::getR() const {return m_gridR->getR();}
const std::vector<double>& SpidaRCVT::getSR() const {return m_gridR->getSR();}
const std::vector<double>& SpidaRCVT::getT() const {return m_gridT->getT();}
const std::vector<double>& SpidaRCVT::getST() const {return m_gridT->getST();}
void SpidaRCVT::RT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->RT_To_SRST(in,out);} 
void SpidaRCVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->SRST_To_RT(in,out);} 
unsigned SpidaRCVT::spectralSize() const {return m_gridR->getNr()*m_gridT->getNt();}
unsigned SpidaRCVT::physicalSize() const {return m_gridR->getNsr()*m_gridT->getNst();}

void SpidaRCVT::mirrorR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const
{
    auto nt = m_gridT->getNt();
    auto nr = m_gridR->getNr();
    for(unsigned j = 0; j < nt; j++){
        for(unsigned i = nr-1; i >= 0; i--)
            out[((nr-1)-i)*nt+j] = in[i*nt+j];
        for(unsigned i = 0; i < nr; i++)
            out[(i+nr)*nt+j] = in[i*nt+j];
    }
}

}
