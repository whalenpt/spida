
#include "spida/ChebInterpX.h"
#include "spida/SpidaCHEBX.h"
#include "spida/grid/chebX.h"
#include "spida/helper/interp.h"
#include "spida/helper/constants.h"
#include <vector>

namespace spida{

  ChebInterpX::ChebInterpX(int ninterp,double minx,double maxx) :
      m_chebx(std::make_unique<SpidaCHEBX>(ChebRootGridX{ninterp,minx,maxx})),
      m_ycheb(ninterp),
      m_dycheb(ninterp) { }

  ChebInterpX::~ChebInterpX() {}

  double ChebInterpX::dXInterp(const std::vector<double>& xin,\
          const std::vector<double>& yin,double xout,int n) 
  {
      SplineInterp spline_to_cheb(xin,yin);
      spline_to_cheb.eval(m_chebx->getGridX().getX(),m_ycheb);
      m_chebx->dX(m_ycheb,m_dycheb,n);
      SplineInterp spline_from_cheb(m_chebx->getGridX().getX(),m_dycheb);
      return spline_from_cheb.eval(xout);
  }

  void ChebInterpX::dXInterp(const std::vector<double>& xin,const std::vector<double>& yin,
          const std::vector<double>& xout,std::vector<double>& dyout,int n)
  {
      SplineInterp spline_to_cheb(xin,yin);
      spline_to_cheb.eval(m_chebx->getGridX().getX(),m_ycheb);
      m_chebx->dX(m_ycheb,m_dycheb,n);
      SplineInterp spline_from_cheb(m_chebx->getGridX().getX(),m_dycheb);
      return spline_from_cheb.eval(xout,dyout);
  }
    
}








