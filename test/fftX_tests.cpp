
#include <gtest/gtest.h>
#include <spida/shape/shapeX.h>
#include <spida/transform/fftCVX.h>
#include <spida/grid/uniformCVX.h>
#include <spida/SpidaCVX.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <random>

// FFTCVX defined such that F{f(x)} = \integral_{-\inf}^{\inf}f(x)exp(-i*kx*x) dx
// Test that forward fft followed by inverse fft yields identity
TEST(FFTCVX_TEST,INVERSES)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);
    for(unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::FFTCVX tr(spida::UniformGridCVX{N,-1,1});
    tr.X_To_SX(in,out);
    tr.SX_To_X(out,expect);

    dataio.writeFile("fft_check.dat",in,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

// F{exp(-a*x^2)} = sqrt(pi/a)*exp(-kx^2/(4*a)) 
TEST(FFTCVX_TEST,GAUSS)
{
	unsigned N = 64;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridCVX grid(N,xmin,xmax);
    const std::vector<double> x = grid.getX();
    const std::vector<double> kx = grid.getSX();
    
    double a = 2.0;
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(-a*pow(x[i],2));
    for(auto i = 0; i < kx.size(); i++)
        expect[i] = sqrt(PI/a)*exp(-pow(kx[i],2)/(4.0*a));

    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);

    dataio.writeFile("fft_gauss.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

// F{sech(a*x)} = (pi/a)*sech(pi*kx/(2*a)) 
TEST(FFTCVX_TEST,SECH)
{
	unsigned N = 64;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridCVX grid(N,xmin,xmax);
    const std::vector<double> x = grid.getX();
    const std::vector<double> kx = grid.getSX();
    
    double a = 2.0;
    for(auto i = 0; i < x.size(); i++)
        in[i] = 1.0/cosh(a*x[i]);
    for(auto i = 0; i < kx.size(); i++)
        expect[i] = (PI/a)/cosh(PI*kx[i]/(2.0*a));

    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);

    dataio.writeFile("fft_sech.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

// F{cos(ax)} = PI*(delta(kx-a) + delta(kx+a))
TEST(FFTCVX_TEST,COS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    spida::UniformGridCVX grid(N,0,2.0*PI);
    const std::vector<double> x = grid.getX();
    const std::vector<double> kx = grid.getSX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = cos(8*x[i]);

    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);
    dataio.writeFile("fft_cos_check_xkx.dat",x,kx);
    dataio.writeFile("fft_cos_check.dat",out);
	EXPECT_DOUBLE_EQ(out[8].real(),PI);
}


TEST(FFTCVX_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");

    using spida::dcmplx;
    using spida::PI;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVX grid(N,0,4*PI);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = sin(x[i]);
    for(auto i = 0; i < x.size(); i++)
        expect[i] = cos(x[i]);

    spida::SpidaCVX spidaX{grid};
    spidaX.dX(in,out);
    dataio.writeFile("fft_der_sin.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

TEST(FFTCVX_TEST,DERIVATIVE_GAUSS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");

    using spida::dcmplx;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVX grid(N,-6,6);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(-pow(x[i],2));
    for(auto i = 0; i < x.size(); i++)
        expect[i] = -2.0*x[i]*exp(-pow(x[i],2));

    spida::SpidaCVX spidaX{grid};
    spidaX.dX(in,out);
    dataio.writeFile("fft_der_gauss.dat",expect,out);
    dataio.writeFile("fft_x_sx.dat",grid.getX(),grid.getSX());
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}





