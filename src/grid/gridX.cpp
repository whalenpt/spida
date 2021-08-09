
#include "spida/grid/gridX.h"
#include <vector>
#include <stdexcept>

// assumes vec is sorted in an increasing order
int indexFromVal(double val,std::vector<double> vec)
{
    if(val < vec[0]){
        std::string msg = "getIndexFromVal(double,vector<double>) error: the val of "\
          + std::to_string(val) + " is less than the minimum xval of " + std::to_string(vec[0]) + ".";
        throw std::invalid_argument(msg);
    } else if(val > vec.back()){
        std::string msg = "getIndexFromVal(double,vector<double>) error: the val of "\
            + std::to_string(val) + " is greater than the maximum xval of " + std::to_string(vec.back()) + ".";
        throw std::invalid_argument(msg);
    }
    auto it = std::upper_bound(std::cbegin(vec),std::cend(vec),val,std::greater<double>());
    if(it == std::begin(vec))
        return 0;
    double diff1 = abs(val - *it);
    double diff2 = abs(val - *(it-1));
    if(diff2 < diff1)
        --it;
    return std::distance(std::cbegin(vec),it);
}

int GridX::getIndexX(double xval) const
{
  if(xval < getMinX()){
      std::string msg = "spida::ChebGridX::getIndexX(double val) error: the val of "\
          + std::to_string(xval) + " is less than the minimum xval of " + std::to_string(getMinX()) + ".";
      throw std::invalid_argument(msg);
  }
  if(xval > getMaxX()){
      std::string msg = "spida::ChebGridX::getIndexX(double val) error: the val of "\
          + std::to_string(xval) + " is greater than the maximum xval of " + std::to_string(getMaxX()) + ".";
      throw std::invalid_argument(msg);
  }
  return indexFromVal(xval,getX());
}








