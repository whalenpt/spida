/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Townsend Whalen   
 *    Email: whalenpt@gmail.com
 *    Date: 09/01/2021
 *    Description: Example of indexing function
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <pwutils/pwindexing.hpp>
#include <iostream>
#include <stdexcept>

//------------------------------------------------------------------------------

template<typename T>
void printNearest(std::vector<T> sorted_vec,T val,unsigned int indx)
{
    std::string array_string = "{";
    for(const auto& item  : sorted_vec)
        array_string += std::to_string(item) + ",";
    array_string.back() = '}';
    std::cout << "The nearest index of " << std::to_string(val) << \
        " in the array " << array_string << " is "\
        + std::to_string(indx) << std::endl;
}

int main()
{
    int int_val = 4;
    std::vector<int> int_vec{0,2,4,6,8,10,12};
    unsigned int int_indx = pw::nearestIndex(int_vec,int_val);
    printNearest(int_vec,int_val,int_indx);

    double double_val = 3.42;
    std::vector<double> double_vec{0.0,1.0,2.0,3.0,4.0,5.0};
    unsigned int double_indx = pw::nearestIndex(double_vec,double_val);
    printNearest(double_vec,double_val,double_indx);

    double double_val2 = 5.0;
    unsigned int double_indx2 = pw::nearestIndex(double_vec,double_val2);
    printNearest(double_vec,double_val2,double_indx2);

    double double_val3 = 0.0;
    unsigned int double_indx3 = pw::nearestIndex(double_vec,double_val3);
    printNearest(double_vec,double_val3,double_indx3);

    try{
        int int_val = -3;
        std::vector<int> int_vec{0,2,4,6,8,10,12};
        unsigned int int_indx = pw::nearestIndex(int_vec,int_val);
        printNearest(int_vec,int_val,int_indx);
    } catch(std::out_of_range& e){
        std::cout << e.what() << std::endl;
    }

    return 0;
}







