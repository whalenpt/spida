
#include "spida/grid/gridT.h"
#include "spida/grid/uniformT.h"
#include "spida/grid/uniformX.h"
#include "spida/constants.h"
#include <string>
#include <iostream>
#include <pwutils/pwexcept.h>

namespace spida{

void setUniformT(double minT,double maxT,std::vector<double>& t)
{
    setUniformX(minT,maxT,t);
}
 
void setUniformST(double minST,double maxST,std::vector<double>& st) 
{
    int nst = st.size();
    double dst = (maxST - minST)/(nst-1);
    if(dst <= 0.0){
        std::string msg = "Error in setUniformST: maxST - minST must be positive.";
        throw pw::Exception(msg);
    }
    for(int i = 0; i < nst; i++) st[i] = minST + i*dst; 
}

UniformGridT::UniformGridT(const UniformGridT& grid) : 
    GridT(grid.getNt(),grid.getMinT(),grid.getMaxT()),
    m_t(grid.getNt()), m_st(grid.getNst())
{
    m_minI = grid.getMinI();
    m_maxI = grid.getMaxI();
    const std::vector<double>& t = grid.getT();
    const std::vector<double>& st = grid.getST();
    std::copy(std::cbegin(t),std::cend(t),std::begin(m_t));
    std::copy(std::cbegin(st),std::cend(st),std::begin(m_st));
}


UniformGridT::UniformGridT(int nt,double minT,double maxT) : 
    GridT(nt,minT,maxT),
    m_t(nt), m_st(nt/2+1), m_minI(0), m_maxI(nt/2)
{
    setUniformT(minT,maxT,m_t);
    double minST = convertIndxToFreq(m_minI);
    double maxST = convertIndxToFreq(m_maxI);
    setUniformST(minST,maxST,m_st);
    checkGrid();
}


UniformGridT::UniformGridT(int nt,double minT,double maxT,\
        double minST,double maxST) : 
    GridT(nt,minT,maxT),
    m_t(nt)
{
    setUniformT(minT,maxT,m_t);
    m_minI = convertFreqToIndx(minST);
    m_maxI = convertFreqToIndx(maxST);
    m_st.resize(m_maxI-m_minI+1,0.0);
    setUniformST(minST,maxST,m_st);
    checkGrid();
}

void UniformGridT::checkMinFreq(double minST)
{
  if(minST < 0.0) {
      std::string msg = "Not a valid minimum frequency of " + std::to_string(getMinST()) +\
          ". Frequencies must be positive.";
      throw pw::Exception(msg);
  }
}

void UniformGridT::checkMaxFreq(double maxST)
{
    // Temporal grid restriction on the maximum frequency (The numerical grid
    // can resolve frequencies up to a certain maximum cutoff based on the # of points in
    // the grid)
    double gridMaxFreq = (2.0*PI/getLT())*(getNt()/2);
    if(maxST > gridMaxFreq)
    {
        auto str1 = std::to_string(gridMaxFreq);
        auto str2 = std::to_string(static_cast<int>(getLT()*maxST/PI)+1);
        std::string msg = "Temporal grid spacing not small enough to accomodate the maximum"\
                           " grid frequency of " + str1 + ". Increase NT to more than " + str2;
        throw pw::Exception(msg);
    }
}

void UniformGridT::checkRangeFreq(double minST,double maxST)
{
  if(maxST <= minST)
  {
      std::string msg = "Not a valid frequency range: The minimum frequency must be less than"\
             " the maximum frequency.";
      throw pw::Exception(msg);
  }
}

void UniformGridT::checkNumSpectralPoints(int nt,int nST)
{
  // Spectral restriction for real fields (maximum number spectral points
  // is nt/2-1 without any zero padding) 
  if(nST > nt/2+1)
  {
      std::string msg = "Frequency Range Too Large. Number of spectral points, "\
          + std::to_string(nST) + ", must less than or equal to " + std::to_string(nt/2+1) \
          + ". Try increasing Nt or reducing the frequency range.";
      throw pw::Exception(msg);
  }
}

void UniformGridT::checkGrid()
{
    checkMinFreq(getMinST());
    checkMaxFreq(getMaxST());
    checkRangeFreq(getMinST(),getMaxST());
    checkNumSpectralPoints(getNt(),getNst());
}

int UniformGridT::convertFreqToIndx(double omeg) const
{
  if(omeg >= 0.0) {
    double dst = 2.0*PI/getLT();
    return round(omeg/dst);
  } else{ 
      std::string msg = "convertFreqToIndxT only processes positive frequencies ";
      throw pw::Exception(msg);
  }
}

double UniformGridT::convertIndxToFreq(int indx) const
{
  if(indx >= 0) {
    double dst = 2.0*PI/getLT();
    return dst*indx;
  } else{ 
      std::string msg = "convertIndxToFreq index must be greater than or equal to 0";
      throw pw::Exception(msg);

  }
}



}


