#pragma once

#include <complex>
#include <vector>
#include "spida/grid/chebX.h"

namespace spida{

  class ChebTransformX 
  {
    public:
        explicit ChebTransformX(const ChebRootGridX& grid);
        ChebTransformX() = delete;
        ChebTransformX(const ChebTransformX&)=delete;
        ChebTransformX& operator=(const ChebTransformX&)=delete;
        void X_To_SX(const std::vector<double>& in,std::vector<double>& out);
        void SX_To_X(const std::vector<double>& in,std::vector<double>& out);
        ~ChebTransformX() = default;
    private:
        int m_nx;
        std::vector<double> m_uFFT;
  };
}




