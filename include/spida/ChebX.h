
#ifndef SPIDACHEBX_H_
#define SPIDACHEBX_H_

#include <vector>

namespace spida{

class ChebTransformX;
class ChebRootGridX;
class GridX;

class ChebX 
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
        const ChebTransformX& getTransformX();
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


