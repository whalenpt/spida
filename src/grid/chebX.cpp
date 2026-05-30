
#include "spida/grid/chebX.h"

#include "spida/helper/constants.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace spida {

void setChebExtremaX(double a, double b, std::vector<double>& x)
{
    if (a >= b) {
        std::string msg = "Error in setChebExtremaX: interval [a,b] -> " + std::to_string(a) +
                          " must be less than " + std::to_string(b);
        throw std::invalid_argument(msg);
    }
    auto nx = static_cast<unsigned>(x.size());
    double dc = spida::PI / (static_cast<double>(nx) - 1);
    double L = b - a;
    for (unsigned j = 0; j < nx; j++)
        x[nx - j - 1] = (cos(j * dc) * L + b + a) / 2.0;
}

void setChebExtremaSX(std::vector<double>& sx)
{
    double dc = spida::PI / (static_cast<double>(sx.size()) - 1);
    for (size_t j = 0; j < sx.size(); j++)
        sx[j] = static_cast<double>(j) * dc;
}

void setChebRootX(double a, double b, std::vector<double>& x)
{
    if (a >= b) {
        std::string msg = "Error in setChebRootX: interval [a,b] -> " + std::to_string(a) +
                          " must be less than " + std::to_string(b);
        throw std::invalid_argument(msg);
    }

    auto nx = static_cast<unsigned>(x.size());
    double dc = spida::PI / static_cast<double>(nx);
    auto L = b - a;
    for (unsigned j = 0; j < nx; j++)
        x[nx - j - 1] = (cos((j + 0.5) * dc) * L + b + a) / 2.0;
}

void setChebRootSX(std::vector<double>& sx)
{
    auto nsx = static_cast<unsigned>(sx.size());
    double dc = spida::PI / static_cast<double>(nsx);
    for (unsigned j = 0; j < nsx; j++)
        sx[j] = (j + 0.5) * dc;
}

ChebExtremaGridX::ChebExtremaGridX(unsigned nx, double min, double max)
    : ChebGridX(nx, min, max), m_x(nx), m_sx(nx)
{
    setChebExtremaX(min, max, this->m_x);
    setChebExtremaSX(this->m_sx);
}

ChebRootGridX::ChebRootGridX(unsigned nx, double min, double max)
    : ChebGridX(nx, min, max), m_x(nx), m_sx(nx)
{
    setChebRootX(min, max, this->m_x);
    setChebRootSX(this->m_sx);
}

} // namespace spida