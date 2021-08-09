
#ifndef SPIDACHEB_X_H_
#define SPIDACHEB_X_H_

#include <vector>
#include "spida/spidaX.h"

namespace spida{

class ChebTransformX;
class ChebRootGridX;

class ChebX : public SpidaX
{
    public:
        ChebX(const ChebRootGridX& grid);
        ChebX() = delete;
        ~ChebX();
        void X_To_SX(const std::vector<double>& in,std::vector<double>& out); 
        void SX_To_X(const std::vector<double>& in,std::vector<double>& out);
        void dX(const std::vector<double>& in,std::vector<double>& out,int n); 
        void dX(const std::vector<double>& in,std::vector<double>& out); 
        const GridX& getGridX(); 
        const TransformX& getTransformX();
    private:
        ChebRootGridX* m_gr;
        ChebTransformX* m_tr;
        std::vector<double> m_sp;
        std::vector<double> m_dsp;
        void dSX(const std::vector<double>& in,std::vector<double>& out,int n); 
        void dSX(const std::vector<double>& in,std::vector<double>& out); 
        void dSX(std::vector<double>& a,int n); 
        void dSX(std::vector<double>& a); 
};

}

#endif


