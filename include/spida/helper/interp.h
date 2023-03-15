
#include <vector>

namespace spida{

// tridisolve solves the tridiagonal linear system Mx=d for x of size N, where a is the
// subdiagonal of size N-1, b is the diagonal of size N, c is the superdiagonal of size N-1,
// and d is the right hand side of size N-1
void tridisolve(const std::vector<double>& a,const std::vector<double>&b,\
        const std::vector<double>&c,const std::vector<double>&d,std::vector<double>&x);
void checkXInterp(const std::vector<double>& xinterp,const std::vector<double>& xdata);
void checkXInterp(double xinterp,const std::vector<double>& xdata);
void checkData(const std::vector<double>& xdata,const std::vector<double>& ydata);

// Linear interpolation class
class LinearInterp{
    public:
      LinearInterp(const std::vector<double>& xvec,const std::vector<double>& yvec);
      std::vector<double> eval(const std::vector<double>& xinterp) const;
      void eval(const std::vector<double>& xinterp,std::vector<double>& yinterp) const;
      double eval(double xinterp) const;

    private:
      std::vector<double> m_xvec;
      std::vector<double> m_yvec;
      void computeInterpY(const std::vector<double>& xinterp,std::vector<double>& yinterp) const;
};

// Cubic spline interpolation class (continuous second derivatives)
class SplineInterp{
    public:
      SplineInterp(const std::vector<double>& xvec,const std::vector<double>& yvec);
      std::vector<double> eval(const std::vector<double>& xinterp) const;
      void eval(const std::vector<double>& xinterp,std::vector<double>& yinterp) const;
      double eval(double xinterp) const;
    private:
      std::vector<double> m_xvec;
      std::vector<double> m_yvec;
      std::vector<double> m_dk;
      std::vector<double> m_ck;
      std::vector<double> m_bk;

      void initializeCoefficients();
      void computeInterpY(const std::vector<double>& xinterp,std::vector<double>& yinterp) const;
};

}