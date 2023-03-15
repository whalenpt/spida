

#include <stdexcept>
#include <string>
#include <pwutils/pwindexing.hpp>
#include "spida/grid/uniformRVT.h"

namespace spida{

std::vector<double> UniformGridRVT::constructGridST(unsigned nt,double minT,double maxT) const
{
    std::vector<double> st(nt/2+1);
    double dst = 2.0*PI/(maxT-minT);
    for(unsigned i = 0; i <= nt/2; i++)
        st[i] = i*dst;
    return st;
}

UniformGridRVT::UniformGridRVT(unsigned nt,double minT,double maxT,\
        double minST,double maxST) : 
    UniformGridT(nt,minT,maxT)
{
    verifyFrequencyRange(minST,maxST);
    std::vector<double> full_st = constructGridST(nt,minT,maxT);
    m_minI = pw::nearestIndex(full_st,minST);
    m_maxI = pw::nearestIndex(full_st,maxST);
    m_nst = m_maxI-m_minI+1;
    m_st.resize(m_nst,0.0);
    for(auto i = m_minI; i <= m_maxI; i++)
        m_st[i-m_minI] = full_st[i];
    m_minST = m_st[0];
    m_maxST = m_st.back();
}

UniformGridRVT::UniformGridRVT(unsigned nt,double minT,double maxT) : 
    UniformGridT(nt,minT,maxT),
    m_st(nt/2+1),
    m_nst(nt/2+1),
    m_minST(indxToFreq(0)),
    m_maxST(indxToFreq(nt/2)),
    m_maxI(nt/2)
{
    m_st = constructGridST(nt,minT,maxT);
}

double UniformGridRVT::maxPossibleFreq() const {
    return indxToFreq(getNt()/2);
}

unsigned UniformGridRVT::freqToIndx(double omeg) const
{
    if(omeg < 0.0){
        std::string msg = "freqToIndx only processes positive frequencies ";
        throw std::domain_error(msg);
    } else if (omeg > maxPossibleFreq()){
        std::string msg = "freqToIndx(double omeg) error: The provided frequency of " \
                   + std::to_string(omeg) + " is bigger than the maximum specified frequency of "\
                   + std::to_string(maxPossibleFreq());
        throw std::domain_error(msg);
    }
    double dst = 2.0*PI/getLT();
    return static_cast<unsigned>(round(omeg/dst));
}

double UniformGridRVT::indxToFreq(unsigned indx) const
{
    if(indx > getNt()/2){
        std::string msg = "indxToFreq(indx) error: Specified index of " + std::to_string(indx) +\
                           " which is bigger than nt/2 of " + std::to_string(getNt()/2);
        throw std::domain_error(msg);
    }
    double dst = 2.0*PI/getLT();
    return dst*indx;
}

void UniformGridRVT::verifyFrequencyRange(double minST,double maxST) const
{

    if(minST >= maxST){
        throw std::domain_error("UniformGridRVT minST of " + std::to_string(minST) +\
                " is greater than or equal to the specified maxST of " + std::to_string(maxST));
    }
    else if(minST < 0.0) {
        throw std::domain_error("UniformGridRVT minST of " + std::to_string(minST) +\
                " is less than zero, please increase to a non-negative value.");
    }
    if(maxST > maxPossibleFreq())
    {
        auto str1 = std::to_string(maxPossibleFreq());
        auto str2 = std::to_string(static_cast<int>((getMaxT()-getMinT())*maxST/PI)+1);
        std::string msg = "Temporal grid spacing not small enough to accomodate the maximum"\
                           " grid frequency of " + str1 + ". Increase NT to more than " + str2;
        throw std::domain_error(msg);
    }
}

}