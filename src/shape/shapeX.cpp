
#include <cmath>
#include "spida/shape/shapeX.h"

namespace spida{

double ExpX::compute(double x) const
{
    return ShapeX::amplitude()*exp(ShapeX::width()*(x-ShapeX::offset())+ShapeX::phase());
}

double CosX::compute(double x) const
{
    return ShapeX::amplitude()*cos(ShapeX::width()*(x-ShapeX::offset())+ShapeX::phase());
}

double SinX::compute(double x) const
{
    return ShapeX::amplitude()*sin(ShapeX::width()*(x-ShapeX::offset())+ShapeX::phase());
}



}








