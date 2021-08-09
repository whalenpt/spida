

#include <gtest/gtest.h>
#include <spida/constants.h>
#include <spida/grid/uniformT.h>
#include <spida/shape/shapeT.h>
//#include <spida/transform/fftwT.h>
#include <spida/transform/periodicT.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <algorithm>
#include <numeric>
#include <functional>

TEST(KTRANSFORM_TEST,GAUSST)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    spida::GaussT shape(std::sqrt(I0),tp,omega0);
    spida::UniformGridT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::PeriodicTransformT transform(grid);

    const std::vector<double>& t = grid.getT();
    std::vector<double> y(nt);
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());
    shape.ShapeT::computeReal(t,y);
    transform.T_To_ST(y,ysp);
    transform.ST_To_T(ysp,yinv);

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),33811981379.747528);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}

TEST(KTRANSFORM_TEST,COMPLEX_GAUSST)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    spida::GaussT shape(std::sqrt(I0),tp,omega0);
    spida::UniformGridT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::PeriodicTransformT transform(grid);

    const std::vector<double>& t = grid.getT();
    std::vector<dcmplx> y(nt);
    std::vector<dcmplx> yinv(nt);

    std::vector<dcmplx> ysp(grid.getNst());
    shape.Shape1D::compute(t,y);
    transform.T_To_ST_c(y,ysp);
    transform.ST_To_T_c(ysp,yinv);

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),67623962759.495056);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}







