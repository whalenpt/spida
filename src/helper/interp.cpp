#include <algorithm>
#include <iostream>
#include <pwutils/pwexcept.h>
#include "spida/helper/interp.h"

namespace spida {

    // Check xinterp vector is valid
    void checkXInterp(const std::vector<double>& xinterp,const std::vector<double>& xdata)
    {
        if(xinterp.empty()){
          throw pw::Exception("eval(std::vector<double> xinterp) "
                   "xinterp must be a non-empty vector");
        }
        if(xinterp[0] < xdata[0]){
            std::string str("eval(std::vector<double> xinterp) "
                    "failed: xinterp[0] is less than the class m_xvec[0], increase "
                    "xinterp[0] to at least ");
            throw pw::Exception(str + std::to_string(xdata[0]));
        }
        if(xinterp.back() > xdata.back()){
            std::string str("eval(std::vector<double> xinterp) "
                    "failed: xinterp[size-1] is greater than the class m_xvec[m_size-1], "
                    "decrease xinterp[size-1] to below or equal to ");
            throw pw::Exception(str + std::to_string(xdata.back()));
        }
    }

    void checkXInterp(double xinterp,const std::vector<double>& xdata)
    {
        if(xinterp < xdata[0]){
            std::string str("eval(double xinterp) "
                    "failed: xinterp is less than the class m_xvec[0], increase "
                    "xinterp[0] to at least ");
            throw pw::Exception(str + std::to_string(xdata[0]));
        }
        if(xinterp > xdata.back()){
            std::string str("eval(double xinterp) "
                    "failed: xinterp is greater than the class m_xvec[m_size-1], "
                    "decrease xinterp to below or equal to ");
            throw pw::Exception(str + std::to_string(xdata.back()));
        }
    }


    void checkData(const std::vector<double>& xdata,const std::vector<double>& ydata)
    {
        if(xdata.size() < 2)
          throw pw::Exception("Interp constructor error: xvec size must be at least 2");
        if(xdata.size() != ydata.size())
          throw pw::Exception("Interp constructor error: the xvec and yvec must be of the same vector size");
    }

    LinearInterp::LinearInterp(const std::vector<double>& xvec,const std::vector<double>& yvec) :
        m_xvec(xvec), m_yvec(yvec)
    {
        checkData(xvec,yvec);
    }

    void LinearInterp::computeInterpY(const std::vector<double>& xinterp,std::vector<double>& yinterp) const
    {
        yinterp.clear();
        yinterp.reserve(xinterp.size());
        unsigned j = 0;
        for(auto xval : xinterp){
            while(xval > m_xvec[j])
                j += 1;
            if(xval == m_xvec[j])
                yinterp.push_back(m_yvec[j]);
            else{
                auto xa = m_xvec[j-1];
                auto xb = m_xvec[j];
                auto ya = m_yvec[j-1];
                auto yb = m_yvec[j];
                yinterp.push_back(ya + ((yb-ya)/(xb-xa))*(xval-xa));
            }
        }
    }

    std::vector<double> LinearInterp::eval(const std::vector<double>& xinterp) const{
        // Check xinterp vector is valid
        checkXInterp(xinterp,m_xvec);

        // Proceed with algorithm
        std::vector<double> yinterp;
        computeInterpY(xinterp,yinterp);
        return yinterp;
    }

    void LinearInterp::eval(const std::vector<double>& xinterp,std::vector<double>& yinterp) const{
        // Check xinterp vector is valid
        checkXInterp(xinterp,m_xvec);
        computeInterpY(xinterp,yinterp);
    }

    double LinearInterp::eval(double xinterp) const
    {
        checkXInterp(xinterp,m_xvec);
        for(unsigned i = 0; i < m_xvec.size(); i++){
            if(m_xvec[i] == xinterp)
                return m_yvec[i];
            else if(m_xvec[i] > xinterp){
                double xa = m_xvec[i-1];
                double xb = m_xvec[i];
                double ya = m_yvec[i-1];
                double yb = m_yvec[i];
                double yinterp = ya + ((yb-ya)/(xb-xa))*(xinterp-xa);
                return yinterp;
            }
        }
        throw pw::Exception("xinterp was not found with the xdata range");
    }

   void tridisolve(const std::vector<double>& a,const std::vector<double>&b,\
        const std::vector<double>&c,const std::vector<double>&d,std::vector<double>&x)
   {
       if(a.size() < 2)
          throw pw::Exception("tridisolve(a,b,c,d,x) error: subdiagonal 'a' size must be at least 2");
       if(b.size() < 3)
          throw pw::Exception("tridisolve(a,b,c,d,x) error: diagonal 'b' size must be at least 3");
       if(c.size() < 2)
          throw pw::Exception("tridisolve(a,b,c,d,x) error: superdiagonal 'c' size must be at least 2");
       if(d.size() < 3)
          throw pw::Exception("tridisolve(a,b,c,d,x) error: R.H.S. 'd' size must be at least 3");
       if(a.size() != (b.size()-1)){
           std::string str("tridisolve(a,b,c,d,x) error: subdiagonal 'a' size of " + std::to_string(a.size())\
                   + " must be equal to the size of the diagonal 'b' " + std::to_string(b.size())\
                   + " minus one");
          throw pw::Exception(str);
       }
       if(c.size() != (b.size()-1)){
           std::string str("tridisolve(a,b,c,d,x) error: superdiagonal 'c' size of " + std::to_string(c.size())\
                   + " must be equal to the size of the diagonal 'b' " + std::to_string(b.size())\
                   + " minus one");
          throw pw::Exception(str);
       }
       if(b.size() != d.size()){
          throw pw::Exception("tridisolve(a,b,c,d,x) error: diagonal 'd' size must be "
              "equal to the R.H.S. 'd' size");
       }
       auto sz = b.size();
       x = std::vector<double>(sz,0.0);
       std::vector<double> bhat(sz);
       std::vector<double> dhat(sz);
       bhat[0] = b[0];
       dhat[0] = d[0];
       for(unsigned i = 1; i < sz; i++){
           double mu = a[i-1]/bhat[i-1];
           bhat[i] = b[i] - mu*c[i-1];
           dhat[i] = d[i] - mu*dhat[i-1];
       }
       x[sz-1] = dhat[sz-1]/bhat[sz-1];
       for(int i = static_cast<int>(sz)-2; i>=0; i--)
           x[i] = (dhat[i] - c[i]*x[i+1])/bhat[i];
   }

    SplineInterp::SplineInterp(const std::vector<double>& xvec,const std::vector<double>& yvec) :
        m_xvec(xvec), m_yvec(yvec),
        m_dk(xvec.size()),m_ck(xvec.size()-1),m_bk(xvec.size()-1) 

    {
        checkData(xvec,yvec);
        initializeCoefficients();
    }

    void SplineInterp::initializeCoefficients(){
        auto sz(m_xvec.size());
        std::vector<double> delta(sz-1);
        std::vector<double> h(sz-1);
        std::vector<double> a(sz-1);
        std::vector<double> b(sz);
        std::vector<double> c(sz-1);
        std::vector<double> r(sz);

        for(unsigned i = 0; i < sz-1; i++){
            h[i] = m_xvec[i+1] - m_xvec[i];
            delta[i] = (m_yvec[i+1] - m_yvec[i])/h[i];
        }
        for(unsigned i = 0; i < sz-2; i++)
            a[i] = h[i+1];
        a[sz-2] = h[sz-2] + h[sz-3];
        b[0] = h[1];
        for(unsigned i = 1; i < sz-1; i++)
            b[i] = 2*(h[i-1] + h[i]);
        b[sz-1] = h[sz-3];
        c[0] = h[0] + h[1];
        for(unsigned i = 1; i < sz-1; i++)
            c[i] = h[i-1];
        r[0] = ((h[0]+2*c[0])*h[1]*delta[0] + h[0]*h[0]*delta[1])/c[0];
        for(unsigned i = 1; i < sz-1; i++)
            r[i] = 3*(h[i-1]*delta[i] + h[i]*delta[i-1]);
        r[sz-1] = (h[sz-2]*h[sz-2]*delta[sz-3]+(2*a[sz-2]+h[sz-2])*h[sz-3]*delta[sz-2])/a[sz-2];

        tridisolve(a,b,c,r,m_dk);
        for(unsigned i = 0; i < sz-1; i++){
            m_ck[i] = (3*delta[i] - 2*m_dk[i] - m_dk[i+1])/h[i];
            m_bk[i] = (m_dk[i] - 2*delta[i] + m_dk[i+1])/(h[i]*h[i]);
        }
    }

    void SplineInterp::computeInterpY(const std::vector<double>& xinterp,std::vector<double>& yinterp) const
    {
        yinterp.clear();
        yinterp.reserve(xinterp.size());

        unsigned j = 0;
        for(auto xval : xinterp){
            while(xval > m_xvec[j])
                j += 1;
            if(xval == m_xvec[j])
                yinterp.push_back(m_yvec[j]);
            else{
                double s = xval - m_xvec[j-1];
                yinterp.push_back(m_yvec[j-1] + s*m_dk[j-1]+s*s*m_ck[j-1] +s*s*s*m_bk[j-1]);
            }
        }
    }


    std::vector<double> SplineInterp::eval(const std::vector<double>& xinterp) const{
        // Check xinterp vector is valid
        checkXInterp(xinterp,m_xvec);
        std::vector<double> yinterp;
        computeInterpY(xinterp,yinterp);
        return yinterp;
    }

    void SplineInterp::eval(const std::vector<double>& xinterp,std::vector<double>& yinterp) const{
        // Check xinterp vector is valid
        checkXInterp(xinterp,m_xvec);
        computeInterpY(xinterp,yinterp);
    }

    double SplineInterp::eval(double xinterp) const
    {
        checkXInterp(xinterp,m_xvec);
        for(unsigned i = 0; i < m_xvec.size(); i++){
            if(m_xvec[i] == xinterp)
                return m_yvec[i];
            else if(m_xvec[i] > xinterp){
                double s = xinterp - m_xvec[i-1];
                double yinterp = m_yvec[i-1] + s*m_dk[i-1]+s*s*m_ck[i-1] +s*s*s*m_bk[i-1];
                return yinterp;
            }
        }
        throw pw::Exception("xinterp was not found with the xdata range");
    }
}
