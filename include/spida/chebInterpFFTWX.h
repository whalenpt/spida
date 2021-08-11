
#ifndef SPIDAINTERPFFTW_X_H_
#define SPIDAINTERPFFTW_X_H_

#include <vector>

namespace spida{
class ChebExtremaGridX;
class ChebFFTWX;

class ChebInterpFFTWX  
{
    public:
        ChebInterpFFTWX(const ChebExtremaGridX& grid);
        ChebInterpFFTWX() = delete;
        ~ChebInterpFFTWX();
        void dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                std::vector<double>& dyout,int n = 1);
        void dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                const std::vector<double>& xout,std::vector<double>& dyout,int n = 1);
        double dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
                double xout,int n = 1);

//        void iX(const std::vector<double>& x,const std::vector<double>& y,std::vector<double>& yout,int n = 1);
//        double iX(const std::vector<double>& x,const std::vector<double>& y,double xmin,double xmax,int n=1);
    private:
        ChebExtremaGridX* m_gr;
        ChebFFTWX* m_chebfftwx;
        std::vector<double> m_ycheb;
        std::vector<double> m_dycheb;
        void interpVec(const std::vector<double>& xin,const std::vector<double>& yin,
                const std::vector<double>& xout,std::vector<double>& yout);
        double interpVal(const std::vector<double>& xin,const std::vector<double>& yin,
                double xout);

};



}

#endif


