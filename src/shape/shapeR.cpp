
#include <cmath>
#include "spida/shape/shapeR.h"

namespace spida{

void compute(const GridR& grid,const ShapeR& shape,std::vector<double>& out)
{
    const std::vector<double>& r = grid.getR();
    for(auto i = 0; i < r.size(); i++)
        out[i] = shape.compute(r[i]);
}


void ShapeR::compute(const GridR& grid,std::vector<double>& y) const
{
    const std::vector<double>& r = grid.getR();
    Shape1D<double,double>::compute(r,y);
}

double GaussR::compute(double r) const
{
    return GaussR::amplitude()*exp(-pow((r-ShapeR::offset())/ShapeR::width(),2));
}

}








