

#include <gtest/gtest.h>
#include <spida/constants.h>
#include <spida/grid/uniformT.h>
#include <spida/grid/besselR.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/periodicT.h>
#include <spida/transform/hankelR.h>
#include <spida/transform/hankelperiodicRT.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <random>

TEST(PERIODICT_TRANSFORM_TEST,GAUSST)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;

    spida::UniformGridT grid(nt,minT,maxT,1.10803e14,1.448963e16);

    std::vector<double> y(nt);
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());

    spida::PeriodicTransformT transform(grid);

    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);
    shape.shapeReal(y);

    transform.T_To_ST(y,ysp);
    transform.ST_To_T(ysp,yinv);

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),33811981379.747528);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}

TEST(PERIODICT_TRANSFORM_TEST,GAUSST_POINTERS)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;

    spida::UniformGridT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::PeriodicTransformT transform(grid);

    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);

    std::vector<double> y(nt);
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());

    shape.ShapeT::shapeReal(y);
    transform.T_To_ST(y.data(),ysp.data());
    transform.ST_To_T(ysp.data(),yinv.data());

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),33811981379.747528);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}


TEST(PERIODICT_TRANSFORM_TEST,COMPLEX_GAUSST)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;

    spida::UniformGridT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::PeriodicTransformT transform(grid);

    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);

    std::vector<dcmplx> y(nt);
    std::vector<dcmplx> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());
    shape.shape(y);

    transform.T_To_ST_c(y,ysp);
    transform.ST_To_T_c(ysp,yinv);

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),67623962759.495056);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}

TEST(HANKEL_TRANSFORM_TEST,GAUSS)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);

    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    double a = 5.0;
    for(auto i = 0; i < N; i++)
        in[i] = exp(-pow(a*r[i],2));
    for(auto i = 0; i < N; i++){
        double beta = 1.0/(2.0*pow(a,2));
        exact[i] = beta*exp(-pow(kr[i],2)/(4.0*pow(a,2)));
    }
    transform.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);
}

TEST(HANKEL_TRANSFORM_TEST,EXP_OVER_R)
{
    int N = 64;
    double rmax = 6.0;
    double a = 1.0;

    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);
    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    for(auto i = 0; i < N; i++)
        in[i] = exp(-a*r[i])/r[i];
    for(auto i = 0; i < N; i++){
        exact[i] = 1.0/sqrt(pow(a,2)+pow(kr[i],2));
    }
    transform.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(exact,out),0.2);
}

TEST(HANKEL_TRANSFORM_TEST,SINC_TEST)
{
    int N = 256;
    double rmax = 32.0;
    double a = 5.0;

    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);
    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    for(auto i = 0; i < N; i++)
        in[i] = r[i] < 1e-4 ? 1.0-pow(a*r[i],2) : sin(a*r[i])/(a*r[i]);
    for(auto i = 0; i < N; i++)
        exact[i] = kr[i] < a ? 1.0/(pow(a,2)*sqrt(1.0-pow(kr[i]/a,2))) : 0.0;

    transform.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(exact,out),0.3);
}



TEST(HANKEL_TRANSFORM_TEST,INVERSES)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR tr(grid);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);
    for(unsigned int i = 0; i < N; i++)
        in[i] = distribution(generator);
    tr.R_To_SR(in,out);
    tr.SR_To_R(out,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

TEST(HANKEL_TRANSFORM_TEST,ORTHOGONALITY)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR tr(grid);

    const std::vector<double>& Ymk = tr.getYmk();
    double zero_sum = 0.0;
    for(auto i = 0; i < N; i++){
        for(auto m = 0; m < N; m++){
            double sum = 0.0;
            for(auto k = 0; k < N; k++){
                sum += Ymk[i*N+k]*Ymk[k*N+m];
            }
            if(i == m)
                EXPECT_NEAR(sum,1,1e-6);
            else
                zero_sum += sum;
        }
    }
    EXPECT_NEAR(zero_sum,0,1e-6);
}

#include <fstream>
TEST(HANKELPERIODICRT_TRANSFORM_TEST,GAUSSTGAUSSR)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 2.0e-15;
    double omega0 = 2.7091e15;
    double minT = -240e-15;
    double maxT = 240e-15;


    spida::UniformGridT gridT(nt,minT,maxT,5.10803e14,1.448963e16);
    spida::GaussT shapeT(gridT,std::sqrt(I0),tp);
    shapeT.setFastPhase(omega0);

    std::vector<double> u0t(nt);
    shapeT.shapeReal(u0t);

    int nr = 100;
    double w0 = 20.0e-6;

    spida::BesselRootGridR gridR(nr,6*w0);
    spida::GaussR shapeR(gridR,1.0,w0);

    std::vector<double> u0r(nr);
    shapeR.shape(u0r);

    std::vector<double> u(nr*nt);
    int nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<double> uinv(nr*nt);

    spida::HankelPeriodicTransformRT transform(gridR,gridT);
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,uinv);

    EXPECT_LT(pw::relative_error(u,uinv),1e-6);
}







