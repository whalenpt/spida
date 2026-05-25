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
        explicit ChebExtremaGridX(unsigned nx,double min=-1.0,double max=1.0);
        ChebExtremaGridX() = delete;
        ~ChebExtremaGridX() override = default;
        const std::vector<double>& getX() const final {return m_x;}
        const std::vector<double>& getSX() const final {return m_sx;}
        double getMinSX() const final {return 0.0;}
        double getMaxSX() const final {return spida::PI;}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

class ChebRootGridX : public ChebGridX
{
    public:
        explicit ChebRootGridX(unsigned nx,double min=-1.0,double max=1.0);
        ChebRootGridX() = delete;
        ~ChebRootGridX() override = default;
        const std::vector<double>& getX() const final {return m_x;}
        const std::vector<double>& getSX() const final {return m_sx;}
        double getMinSX() const final {return 0.0;}
        double getMaxSX() const final {return spida::PI;}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

}