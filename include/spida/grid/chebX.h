#pragma once 

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
        using GridX::GridX;
        ~ChebGridX() override = default; 
};

class ChebExtremaGridX : public ChebGridX
{
    public:
        // Dont allow implicit conversion from ChebTransformFFTWX(int nx) to ChebTransform(ChebExtremaGridX)
        explicit ChebExtremaGridX(int nx,double min=-1.0,double max=1.0);
        ChebExtremaGridX() = delete;
        ~ChebExtremaGridX() final = default;
        const std::vector<double>& getX() const final {return m_x;}
        const std::vector<double>& getSX() const final {return m_sx;}
        double getL() const {return getMaxX()-getMinX();}
        double getMinSX() const final {return 0.0;}
        double getMaxSX() const final {return spida::PI;}
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
        ~ChebRootGridX() final = default;
        const std::vector<double>& getX() const final {return m_x;}
        const std::vector<double>& getSX() const final {return m_sx;}
        double getL() const {return getMaxX()-getMinX();}
        double getMinSX() const final {return 0.0;}
        double getMaxSX() const final {return spida::PI;}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

}