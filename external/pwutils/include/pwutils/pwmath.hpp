
#ifndef PWMATH_HPP_ 
#define PWMATH_HPP_

#include <string>
#include <vector>
#include <cmath> //sqrt
#include <numeric> //accumulate
#include <algorithm> //max_element, min_element
#include <complex>

namespace pw{

    template<typename T>
    auto norm(const std::vector<std::complex<T>>& v) {
        T squared_sum = std::accumulate(std::begin(v),std::end(v),0.0,\
                [](auto sum,auto& val) {return sum + std::norm(val);});
        return sqrt(squared_sum);
    }

    template<typename T>
    auto max_element(const std::vector<std::complex<T>>& v) {
        auto dcmplx_lt = [](auto a,auto b) {return abs(a) < abs(b);};
        return std::max_element(std::begin(v),std::end(v),dcmplx_lt);
    }

    template<typename T>
    auto max(const std::vector<std::complex<T>>& v) {
        auto itmax = pw::max_element<T>(v);
        return *itmax;
    }

    template<typename T>
    auto argmax(const std::vector<std::complex<T>>& v) {
        auto itmax = pw::max_element<T>(v);
        return std::distance(std::begin(v),itmax);
    }

    template<typename T>
    auto min_element(const std::vector<std::complex<T>>& v) {
        auto dcmplx_lt = [](auto a,auto b) {return abs(a) < abs(b);};
        return std::min_element(std::begin(v),std::end(v),dcmplx_lt);
    }
    template<typename T>
    auto min(const std::vector<std::complex<T>>& v) {
        auto itmin = pw::min_element<T>(v);
        return *itmin;
    }
    template<typename T>
    auto argmin(const std::vector<std::complex<T>>& v) {
        auto itmin = pw::min_element<T>(v);
        return std::distance(std::begin(v),itmin);
    }

    template<typename T>
    auto relative_error(const std::vector<std::complex<T>>& v1,\
            const std::vector<std::complex<T>>& v2){
        assert(v1.size() == v2.size());
        auto v1_norm = pw::norm<T>(v1);
        T sum = 0;
        for(auto i = 0; i < v1.size(); i++)
            sum += std::norm(v1[i]-v2[i]); 
        return sqrt(sum)/v1_norm;
    }

    template<typename T>
    auto norm(const std::vector<T>& v) {
        T squared_sum = std::accumulate(std::begin(v),std::end(v),0.0,\
                [](auto sum,auto& val) {return sum + pow(val,2);});
        return sqrt(squared_sum);
    }

    template<typename T>
    auto max(const std::vector<T>& v) {
        return *std::max_element(std::begin(v),std::end(v));
    }

    template<typename T>
    auto argmax(const std::vector<T>& v) {
        auto itmax = std::max_element<T>(std::begin(v),std::end(v));
        return std::distance(std::begin(v),itmax);
    }

    template<typename T>
    auto min(const std::vector<T>& v) {
        return *std::min_element(std::begin(v),std::end(v));
    }

    template<typename T>
    auto argmin(const std::vector<T>& v) {
        auto itmin = std::min_element<T>(std::begin(v),std::end(v));
        return std::distance(std::begin(v),itmin);
    }

    template<typename T>
    auto relative_error(const std::vector<T>& v1,const std::vector<T>& v2){
        assert(v1.size() == v2.size());
        auto v1_norm = pw::norm<T>(v1);
        T sum = 0;
        for(auto i = 0; i < v1.size(); i++)
            sum += pow(v1[i]-v2[i],2); 
        return sqrt(sum)/v1_norm;
    }


    int factorial(int);
    bool isInteger(const std::string& s);
    bool rowIsIntegers(const std::vector<std::string>& row);
    bool lineIsIntegers(const std::string& s);
    bool isDouble(const std::string& s);
    bool rowIsDoubles(const std::vector<std::string>& row);
    bool lineIsDoubles(const std::string& s);
    int intceil(int x,int y);
    unsigned int intceil(unsigned int x,unsigned int y);
}

#endif


