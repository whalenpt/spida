 
#ifndef SPIDA_GRID_CHEBX_H_
#define SPIDA_GRID_CHEBX_H_

#include <vector>
#include "spida/grid/gridX.h"
#include "spida/helper/constants.h"

namespace spida{

void setChebExtremaX(double minX,double maxX,std::vector<double>& x);
void setChebExtremaSX(std::vector<double>& sx);
void setChebRootX(double minX,double maxX,std::vector<double>& x);
void setChebRootSX(std::vector<double>& sx);


class ChebGridX : public GridX
{
    public:
        ChebGridX(int nx,double min,double max) :
            GridX(nx,min,max) {}
        ~ChebGridX() {}; 
        virtual const std::vector<double>& getX() const = 0;
        virtual const std::vector<double>& getSX() const = 0;
        virtual double getMinSX() const = 0;
        virtual double getMaxSX() const = 0;
};

class ChebExtremaGridX : public ChebGridX
{
    public:
        // Dont allow implicit conversion from ChebTransformFFTWX(int nx) to ChebTransform(ChebExtremaGridX)
        explicit ChebExtremaGridX(int nx,double min=-1.0,double max=1.0);
        ChebExtremaGridX() = delete;
        ~ChebExtremaGridX() {}; 
        const std::vector<double>& getX() const {return m_x;}
        const std::vector<double>& getSX() const {return m_sx;}
        double getL() const {return getMaxX()-getMinX();}
        double getMinSX() const {return 0.0;}
        double getMaxSX() const {return spida::PI;}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

class ChebRootGridX : public ChebGridX
{
    public:
        // Dont allow implicit conversion from ChebTransformX(int nx) to ChebTransform(ChebRootGridX)
        explicit ChebRootGridX(int nx,double min=-1.0,double max=1.0);
        ChebRootGridX() = delete;

        ~ChebRootGridX() {}; 
        const std::vector<double>& getX() const {return m_x;}
        const std::vector<double>& getSX() const {return m_sx;}
        double getL() const {return getMaxX()-getMinX();}
        double getMinSX() const {return 0.0;}
        double getMaxSX() const {return spida::PI;}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};


}

#endif


