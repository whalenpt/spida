
#ifndef SPIDA_GRID_BESSELR_H_
#define SPIDA_GRID_BESSELR_H_

#include <vector>
#include "spida/grid/gridR.h"
#include "spida/constants.h"

namespace spida{

class BesselRootGridR : public GridR
{
    public:
        explicit BesselRootGridR(int nr,double maxr); 
        BesselRootGridR() = delete;
        ~BesselRootGridR() {}; 
        const std::vector<double>& getR() const {return m_r;}
        const std::vector<double>& getSR() const {return m_sr;}
        const std::vector<double>& getBesselRoots() const {return m_roots;}
        double getMaxSR() const {return m_jN/GridR::getMaxR();}
        double getjN() const {return m_jN;}
    private:
        std::vector<double> m_r;
        std::vector<double> m_sr;
        std::vector<double> m_roots;
        double m_jN;
};

}

#endif // SPIDA_GRID_BESSELZEROR_H_


