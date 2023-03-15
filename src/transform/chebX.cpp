#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <nayukidct/FastDctLee.hpp>
#include "spida/grid/chebX.h" 
#include "spida/transform/chebX.h"

namespace spida {

  ChebTransformX::ChebTransformX(const ChebRootGridX& grid) : 
      m_nx(grid.getNx()),
      m_uFFT(grid.getNx()) {}

  void ChebTransformX::X_To_SX(const std::vector<double>& in,std::vector<double>& out) 
  {
        std::reverse_copy(std::begin(in),std::end(in),std::begin(m_uFFT));
        FastDctLee::transform(m_uFFT);
        for(int j = 0; j < m_nx; j++)  
            out[j] = 2.0*m_uFFT[j]/static_cast<double>(m_nx);
  }

  void ChebTransformX::SX_To_X(const std::vector<double>& in,std::vector<double>& out)
  {
        std::copy(std::begin(in),std::end(in),std::begin(m_uFFT));
        FastDctLee::inverseTransform(m_uFFT);
        std::reverse_copy(std::begin(m_uFFT),std::end(m_uFFT),std::begin(out));
  }

}