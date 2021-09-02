
#include "spida/grid/gridX.h"
#include <vector>
#include <stdexcept>
#include <pwutils/pwindexing.hpp>

unsigned int GridX::getIndexX(double xval) const
{
    try{
        unsigned int indx = pw::nearestIndex(getX(),xval);
        return indx;
    } catch(std::exception& e){
      std::string msg = "spida::ChebGridX::getIndexX(double val) error: the val of "\
          + std::to_string(xval) + " was not within the GridX specified range of [" \
          + std::to_string(getMinX()) + "," + std::to_string(getMaxX()) + "]";
      throw std::invalid_argument(msg);
    }
}








