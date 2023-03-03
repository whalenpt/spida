#pragma once

#include <memory>
#include <vector>
#include "spida/chebX.h"

namespace spida{

class ChebInterpX  
{
    public:
        ChebInterpX(int ninterp,double minx,double maxx); 
        ChebInterpX() = delete;
        ~ChebInterpX() = default;
        void dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                const std::vector<double>& xout,std::vector<double>& dyout,int n = 1);
        void dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                std::vector<double>& dyout,int n = 1) {dXInterp(xin,yin,xin,dyout,n);}
        double dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                double xout,int n = 1);
    private:
        std::unique_ptr<SpidaChebX> m_chebx;
        std::vector<double> m_ycheb;
        std::vector<double> m_dycheb;
};

}