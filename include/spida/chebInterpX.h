
#ifndef CHEBINTERPX_H_
#define CHEBINTERPX_H_

#include <vector>

namespace spida{
class ChebX;

class ChebInterpX  
{
    public:
        ChebInterpX(int ninterp,int minx,int maxx); 
        ChebInterpX() = delete;
        ~ChebInterpX();
        void dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                const std::vector<double>& xout,std::vector<double>& dyout,int n = 1);
        void dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                std::vector<double>& dyout,int n = 1) {dXInterp(xin,yin,xin,dyout,n);}
        double dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                double xout,int n = 1);
//        void iX(const std::vector<double>& x,const std::vector<double>& y,std::vector<double>& yout,int n = 1);
//        double iX(const std::vector<double>& x,const std::vector<double>& y,double xmin,double xmax,int n=1);
    private:
        ChebX* m_chebx;
        std::vector<double> m_ycheb;
        std::vector<double> m_dycheb;
};



}

#endif


