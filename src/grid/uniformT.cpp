
#include "spida/grid/gridT.h"
#include "spida/grid/uniformT.h"
#include "spida/helper/constants.h"
#include <string>
#include <iostream>
#include <stdexcept>
#include <pwutils/pwindexing.hpp>

namespace spida{

std::vector<double> buildUniformT(unsigned nt,double minT,double maxT)
{
    if(minT >= maxT){
        std::string msg = "Error in buildUniformT: minT must be less than maxT.";
        throw std::invalid_argument(msg);
    }
    std::vector<double> t(nt);
    double dt = (maxT - minT)/static_cast<double>(nt-1);
    for(auto i = 0; i < nt; i++) t[i] = minT + i*dt; 
    return t;
}

std::vector<double> buildUniformST(unsigned nt,double minT,double maxT)
{
    if(minT >= maxT){
        std::string msg = "Error in buildUniformST: minT must be less than maxT.";
        throw std::invalid_argument(msg);
    }
    std::vector<double> full_st(nt/2+1);
    double dst = 2.0*PI/(maxT-minT);
    for(auto i = 0; i <= nt/2; i++)
        full_st[i] = i*dst;
    return full_st;
}

UniformGridT::UniformGridT(const UniformGridT& grid) : 
    GridT(grid.getNt(),grid.getMinT(),grid.getMaxT()),
    m_nst(grid.getNst()),
    m_minST(grid.getMinST()),
    m_maxST(grid.getMaxST()),
    m_minI(grid.getMinI()),
    m_maxI(grid.getMaxI()),
    m_t(grid.getNt()), 
    m_st(grid.getNst())
{
    const std::vector<double>& t = grid.getT();
    const std::vector<double>& st = grid.getST();
    std::copy(std::cbegin(t),std::cend(t),std::begin(m_t));
    std::copy(std::cbegin(st),std::cend(st),std::begin(m_st));
}

UniformGridT::UniformGridT(unsigned int nt,double minT,double maxT) : 
    GridT(nt,minT,maxT),
    m_nst(nt/2+1),
    m_minST(indxToFreq(0)),
    m_maxST(indxToFreq(nt/2)),
    m_minI(0),
    m_maxI(nt/2)
{
    m_t = buildUniformT(nt,minT,maxT);
    m_st = buildUniformST(nt,minT,maxT);
}


UniformGridT::UniformGridT(unsigned int nt,double minT,double maxT,\
        double minST,double maxST) : 
    GridT(nt,minT,maxT)
{
    verifyFrequencyRange(minST,maxST);
    m_t = buildUniformT(nt,minT,maxT);
    std::vector<double> full_st = buildUniformST(nt,minT,maxT);
    m_minI = pw::nearestIndex(full_st,minST);
    m_maxI = pw::nearestIndex(full_st,maxST);
    m_nst = m_maxI-m_minI+1;
    m_st.resize(m_nst,0.0);
    for(auto i = m_minI; i <= m_maxI; i++)
        m_st[i-m_minI] = full_st[i];
    m_minST = m_st[0];
    m_maxST = m_st.back();
}

double UniformGridT::maxPossibleFreq() const {
    return indxToFreq(getNt()/2);
}

unsigned int UniformGridT::freqToIndx(double omeg) const
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
    return round(omeg/dst);
}

double UniformGridT::indxToFreq(unsigned int indx) const
{
    if(indx > getNt()/2){
        std::string msg = "indxToFreq(indx) error: Specified index of " + std::to_string(indx) +\
                           " which is bigger than nt/2 of " + std::to_string(getNt()/2);
        throw std::domain_error(msg);
    }
    double dst = 2.0*PI/getLT();
    return dst*indx;
}

void UniformGridT::verifyFrequencyRange(double minST,double maxST) const
{

    if(minST >= maxST){
        throw std::domain_error("UniformGridT minST of " + std::to_string(minST) +\
                " is greater than or equal to the specified maxST of " + std::to_string(maxST));
    }
    else if(minST < 0.0) {
        throw std::domain_error("UniformGridT minST of " + std::to_string(minST) +\
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


