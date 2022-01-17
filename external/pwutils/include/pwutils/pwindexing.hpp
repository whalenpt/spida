#include <vector>
#include <stdexcept>
#include <string>
#include <cmath>

namespace pw{

// Assume vector array is sorted from lowest to highest. Check first that val is within the vector range."
template<typename T>
unsigned int nearestIndex(std::vector<T>& sorted_vec,T val)
{
    if(val < sorted_vec[0]){
        std::string msg = "nearestIndex(T val,vector<T>& sorted_vec) error: the val of "\
          + std::to_string(val) + " is less than the minimum xval of " + std::to_string(sorted_vec[0]) + ".";
        throw std::out_of_range(msg);
    } else if(val > sorted_vec.back()){
        std::string msg = "nearestIndex(T val,vector<T>& sorted_vec) error: the val of "\
            + std::to_string(val) + " is greater than the maximum xval of " + std::to_string(sorted_vec.back()) + ".";
        throw std::out_of_range(msg);
    }
    auto it = std::lower_bound(std::cbegin(sorted_vec),std::cend(sorted_vec),val);
    if(it == sorted_vec.cbegin())
        return 0;
    unsigned target = std::distance(std::cbegin(sorted_vec),it);

    auto a = *(it - 1);
    auto b = *(it);
    return (std::abs(val - a) < std::abs(val-b) ? (target-1) : target);
}     


}

