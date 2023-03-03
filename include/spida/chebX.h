#pragma once

#include <memory>
#include <vector>
#include "spida/grid/chebX.h"
#include "spida/transform/chebX.h"

namespace spida{

class SpidaChebX 
{
    public:
        explicit SpidaChebX(const ChebRootGridX& grid);
        SpidaChebX() = delete;
        ~SpidaChebX() = default;
        void X_To_SX(const std::vector<double>& in,std::vector<double>& out); 
        void SX_To_X(const std::vector<double>& in,std::vector<double>& out);
        void dX(const std::vector<double>& in,std::vector<double>& out,int n); 
        void dX(const std::vector<double>& in,std::vector<double>& out); 
        const GridX& getGridX(); 
        const ChebTransformX& getTransformX();
    private:
        std::unique_ptr<ChebRootGridX> m_gr;
        std::unique_ptr<ChebTransformX> m_tr;
        std::vector<double> m_sp;
        std::vector<double> m_dsp;
        void dSX(const std::vector<double>& in,std::vector<double>& out,int n); 
        void dSX(const std::vector<double>& in,std::vector<double>& out); 
        void dSX(std::vector<double>& a,int n); 
        void dSX(std::vector<double>& a); 
};

}