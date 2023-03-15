
#include <gtest/gtest.h>
#include <spida/CVX.h>
#include <spida/RVX.h>
#include <spida/shape/shapeX.h>
#include <spida/transform/fftCVX.h>
#include <spida/transform/fftRVX.h>
#include <spida/grid/uniformCVX.h>
#include <spida/grid/uniformRVX.h>
#include <pwutils/pwmath.hpp>
#include <random>

// FFTCVX defined such that F{f(x)} = \integral_{-\inf}^{\inf}f(x)exp(-i*kx*x) dx
// Test that forward fft followed by inverse fft yields identity
TEST(FFTCVX_TEST,INVERSES)
{
	unsigned N = 32;
    using spida::dcmplx;

    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);
    for(unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::FFTCVX tr(spida::UniformGridCVX{N,-1,1});
    tr.X_To_SX(in,out);
    tr.SX_To_X(out,expect);

    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

TEST(FFTCVX_TEST,GAUSS)
{
    // F{exp(-a*x^2)} = sqrt(pi/a)*exp(-kx^2/(4*a)) 
	unsigned N = 64;
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridCVX grid(N,xmin,xmax);
    const auto& x = grid.getX();
    const auto& kx = grid.getSX();
    
    double a = 2.0;
    for(size_t i = 0; i < x.size(); i++)
        in[i] = exp(-a*pow(x[i],2));
    for(size_t i = 0; i < kx.size(); i++)
        expect[i] = sqrt(PI/a)*exp(-pow(kx[i],2)/(4.0*a));

    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);

    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

TEST(FFTCVX_TEST,SECH)
{
    // F{sech(a*x)} = (pi/a)*sech(pi*kx/(2*a)) 
	unsigned N = 64;
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridCVX grid(N,xmin,xmax);
    const auto& x = grid.getX();
    const auto& kx = grid.getSX();
    
    double a = 2.0;
    for(size_t i = 0; i < x.size(); i++)
        in[i] = 1.0/cosh(a*x[i]);
    for(size_t i = 0; i < kx.size(); i++)
        expect[i] = (PI/a)/cosh(PI*kx[i]/(2.0*a));

    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);

    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

TEST(FFTCVX_TEST,COS)
{
    // F{cos(ax)} = PI*(delta(kx-a) + delta(kx+a))
	unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    spida::UniformGridCVX grid(N,0,2.0*PI);
    const auto& x = grid.getX();
    for(size_t i = 0; i < x.size(); i++)
        in[i] = cos(8*x[i]);

    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);
	EXPECT_DOUBLE_EQ(out[8].real(),PI);
}


TEST(FFTCVX_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;

    using spida::dcmplx;
    using spida::PI;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVX grid(N,0,4*PI);
    const auto& x = grid.getX();
    for(size_t i = 0; i < x.size(); i++)
        in[i] = sin(x[i]);
    for(size_t i = 0; i < x.size(); i++)
        expect[i] = cos(x[i]);

    spida::SpidaCVX spidaX{grid};
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

TEST(FFTCVX_TEST,DERIVATIVE_GAUSS)
{
	unsigned N = 32;

    using spida::dcmplx;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVX grid(N,-6,6);
    const auto& x = grid.getX();
    for(size_t i = 0; i < x.size(); i++)
        in[i] = exp(-pow(x[i],2));
    for(size_t i = 0; i < x.size(); i++)
        expect[i] = -2.0*x[i]*exp(-pow(x[i],2));

    spida::SpidaCVX spidaX{grid};
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

// FFTRVX defined such that F{f(x)} = \integral_{-\inf}^{\inf}f(x)exp(-i*kx*x) dx
// Test that forward fft followed by inverse fft yields identity
TEST(FFTRVX_TEST,INVERSES)
{
    using spida::dcmplx;

	unsigned N = 32;
    spida::UniformGridRVX grid{N,-2,2};
    unsigned nsx = grid.getNsx();
    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<double> in(N);
    std::vector<dcmplx> out(nsx);
    std::vector<double> expect(N);
    for(unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::FFTRVX tr(grid);
    tr.X_To_SX(in,out);
    tr.SX_To_X(out,expect);

    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

TEST(FFTRVX_TEST,GAUSS)
{
    // F{exp(-a*x^2)} = sqrt(pi/a)*exp(-kx^2/(4*a)) 
	unsigned N = 64;
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridRVX grid(N,xmin,xmax);
    unsigned nsx =grid.getNsx();
    std::vector<double> in(N,0.0);
    std::vector<dcmplx> out(nsx,0.0);
    std::vector<dcmplx> expect(nsx,0.0);

    const auto& x = grid.getX();
    const auto& kx = grid.getSX();
    
    double a = 2.0;
    for(size_t i = 0; i < x.size(); i++)
        in[i] = exp(-a*pow(x[i],2));
    for(size_t i = 0; i < kx.size(); i++)
        expect[i] = sqrt(PI/a)*exp(-pow(kx[i],2)/(4.0*a));

    spida::FFTRVX tr(grid);
    tr.X_To_SX(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

TEST(FFTRVX_TEST,SECH)
{
    // F{sech(a*x)} = (pi/a)*sech(pi*kx/(2*a)) 
	unsigned N = 64;
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridRVX grid(N,xmin,xmax);
    unsigned nsx =grid.getNsx();
    std::vector<double> in(N,0.0);
    std::vector<dcmplx> out(nsx,0.0);
    std::vector<dcmplx> expect(nsx,0.0);


    const auto& x = grid.getX();
    const auto& kx = grid.getSX();
    
    double a = 2.0;
    for(size_t i = 0; i < x.size(); i++)
        in[i] = 1.0/cosh(a*x[i]);
    for(size_t i = 0; i < kx.size(); i++)
        expect[i] = (PI/a)/cosh(PI*kx[i]/(2.0*a));

    spida::FFTRVX tr(grid);
    tr.X_To_SX(in,out);

    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

TEST(FFTRVX_TEST,COS)
{
    // F{cos(ax)} = PI*(delta(kx-a) + delta(kx+a))
	unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    // Need cos(xmin) = cos(xmax) for periodicity
    double xmin = 0.0;
    double xmax = 2.0*PI;
    spida::UniformGridRVX grid(N,xmin,xmax);
    unsigned nsx = grid.getNsx();
    std::vector<double> in(N);
    std::vector<dcmplx> out(nsx);
    const auto& x = grid.getX();
    for(size_t i = 0; i < x.size(); i++)
        in[i] = cos(8*x[i]);

    spida::FFTRVX tr(grid);
    tr.X_To_SX(in,out);
	EXPECT_DOUBLE_EQ(out[8].real(),PI);
}


TEST(FFTRVX_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    // Need sin(xmin) = sin(xmax) for periodicity
    double xmin = 0.0;
    double xmax = 2.0*PI;
    spida::UniformGridRVX grid(N,xmin,xmax);
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);

    const auto& x = grid.getX();
    for(size_t i = 0; i < x.size(); i++)
        in[i] = sin(x[i]);
    for(size_t i = 0; i < x.size(); i++)
        expect[i] = cos(x[i]);

    spida::SpidaRVX spidaX{grid};
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-4);
}

TEST(FFTRVX_TEST,DERIVATIVE_GAUSS)
{
	unsigned N = 32;
    using spida::dcmplx;

    spida::UniformGridRVX grid(N,-6,6);
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);

    const auto& x = grid.getX();
    for(size_t i = 0; i < x.size(); i++)
        in[i] = exp(-pow(x[i],2));
    for(size_t i = 0; i < x.size(); i++)
        expect[i] = -2.0*x[i]*exp(-pow(x[i],2));

    spida::SpidaRVX spidaX{grid};
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}








