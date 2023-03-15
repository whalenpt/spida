

#include <gtest/gtest.h>
#include <spida/CVT.h>
#include <spida/RVT.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformCVT.h>
#include <spida/shape/shapeT.h>
#include <spida/transform/fftRVT.h>
#include <spida/transform/fftCVT.h>
#include <pwutils/pwmath.hpp>
#include <fstream>
#include <random>


// FFTCVT defined such that F{f(t)} = \integral_{-\inf}^{\inf}f(t)exp(i*omega*t) dt
// Test that forward fft followed by inverse fft yields identity
TEST(FFTCVT_TEST,INVERSES)
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

    spida::FFTCVT tr(spida::UniformGridCVT{N,-1,1});
    tr.T_To_ST(in,out);
    tr.ST_To_T(out,expect);

    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

// a = 1/tp^2
// F{exp(-a*t^2)} = sqrt(pi/a)*exp(-omega^2/(4*a)) 
// F{exp(-(t/tp)^2)}= tp*sqrt(pi)*exp(-tp^2*omega^2/4)
TEST(FFTCVT_TEST,GAUSS)
{
	unsigned N = 64;
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    double xmin = -6;
    double xmax = 6;
    spida::UniformGridCVT grid(N,xmin,xmax);
    const std::vector<double> t = grid.getT();
    const std::vector<double> omega = grid.getST();
    
    double a = 2.0;
    for(size_t i = 0; i < t.size(); i++)
        in[i] = exp(-a*pow(t[i],2));
    for(size_t i = 0; i < omega.size(); i++)
        expect[i] = sqrt(PI/a)*exp(-pow(omega[i],2)/(4.0*a));

    spida::FFTCVT tr(grid);
    tr.T_To_ST(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

TEST(FFTCVT_TEST,COS)
{
    // F{cos(at)} = PI*(delta(omega-a) + delta(omega+a))
	unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    spida::UniformGridCVT grid(N,0.0,2.0*PI);
    const std::vector<double> t = grid.getT();
    for(size_t i = 0; i < t.size(); i++)
        in[i] = cos(8*t[i]);

    spida::FFTCVT tr(grid);
    tr.T_To_ST(in,out);
	EXPECT_DOUBLE_EQ(out[8].real(),PI);
}

TEST(FFTCVT_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;

    using spida::dcmplx;
    using spida::PI;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVT grid(N,0,2*PI);
    const std::vector<double> t = grid.getT();
    for(size_t i = 0; i < t.size(); i++)
        in[i] = sin(t[i]);
    for(size_t i = 0; i < t.size(); i++)
        expect[i] = cos(t[i]);

    spida::SpidaCVT spi{grid};
    spi.dT(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

TEST(FFTCVT_TEST,DERIVATIVE_GAUSS)
{
	unsigned N = 32;

    using spida::dcmplx;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVT grid(N,-6,6);
    const std::vector<double> t = grid.getT();
    for(size_t i = 0; i < t.size(); i++)
        in[i] = exp(-pow(t[i],2));
    for(size_t i = 0; i < t.size(); i++)
        expect[i] = -2.0*t[i]*exp(-pow(t[i],2));

    spida::SpidaCVT spi{grid};
    spi.dT(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

// FFTRVT defined such that F{f(t)} = \integral_{-\inf}^{\inf}f(t)exp(i*omega*t) dt
// Test that forward fft followed by inverse fft yields identity
TEST(FFTRVT_TEST,INVERSES)
{
    using spida::dcmplx;

	unsigned N = 32;
    spida::UniformGridRVT grid{N,-2,2};
    unsigned nst = grid.getNst();

    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<double> in(N);
    std::vector<dcmplx> out(nst);
    std::vector<double> expect(N);
    for(unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::FFTRVT tr(grid);
    tr.T_To_ST(in,out);
    tr.ST_To_T(out,expect);

    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}


// FFT{exp(-(t/tp)^2)exp(-i*omega0*t}= tp*sqrt(pi)*exp(-tp^2*(omega-omega0)^2/4)
TEST(FFTRVT_TEST,GAUSST)
{
    using spida::dcmplx;
    using spida::PI;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;
    double minST = 1.10803e14;
    double maxST = 1.448963e16;

    spida::UniformGridRVT grid(nt,minT,maxT,minST,maxST);
    unsigned nst = grid.getNst();
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(nst);

    spida::FFTRVT transform(grid);
    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);
    auto y = shape.shapeRV();

    transform.T_To_ST(y,ysp);
    transform.ST_To_T(ysp,yinv);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);

    std::vector<dcmplx> ysp_ex(grid.getNst(),0.0);
    const std::vector<double>& omega = grid.getST();
    // y = f(t)*cos(i\omega0t) - > FFT{y} = (FFT{f(\omega - \omega0)}+FFT{f(\omega+\omega0)})/2
    // For real fields, fft taken over positive frequencies: FFT_real{y} = FFT_real{f(\omega-\omega0)}/2  
    for(size_t j = 0; j < grid.getNst(); j++)
        ysp_ex[j] = 0.5*std::sqrt(I0)*tp*sqrt(PI)*exp(-pow(tp,2)*pow(omega[j]-omega0,2)/4.0);

    EXPECT_LT(pw::relative_error(ysp,ysp_ex),1e-5);

    auto ycmplx = shape.shapeCV();
    std::vector<dcmplx> ysp_exCV(grid.getNst(),0.0);
    // y = f(t)*exp(i\omega0t) - > FFT{y} = FFT{f(\omega - \omega0)}
    // For complex fields, multiplication by exp(i\omega0t) in real space is a simple shift in spectral space
    for(size_t j = 0; j < grid.getNst(); j++)
        ysp_exCV[j] = std::sqrt(I0)*tp*sqrt(PI)*exp(-pow(tp,2)*pow(omega[j]-omega0,2)/4.0);

    // Complex valued transform -> phase works fine
    transform.CVT_To_ST(ycmplx,ysp_exCV);
    EXPECT_LT(pw::relative_error(ysp,ysp_ex),1e-6);
}

TEST(FFTRVT_TEST,COS)
{
    // F{cos(at)} = PI*(delta(omega-a) + delta(omega+a))
	unsigned N = 32;
    using spida::dcmplx;
    using spida::PI;

    std::vector<double> in(N);
    spida::UniformGridRVT grid(N,0.0,2.0*PI);
	unsigned nst = grid.getNst();
    std::vector<dcmplx> out(nst);

    const std::vector<double> t = grid.getT();
    for(size_t i = 0; i < t.size(); i++)
        in[i] = cos(8*t[i]);

    spida::FFTRVT tr(grid);
    tr.T_To_ST(in,out);
	EXPECT_DOUBLE_EQ(out[8].real(),PI);
}

TEST(FFTRVT_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;

    using spida::dcmplx;
    using spida::PI;
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);

    // Need sin(tmin) = sin(tmax) for periodicity
    double tmin = 0.0;
    double tmax = 2.0*PI;
    spida::UniformGridRVT grid(N,tmin,tmax);

    const std::vector<double> t = grid.getT();
    for(size_t i = 0; i < t.size(); i++)
        in[i] = sin(t[i]);
    for(size_t i = 0; i < t.size(); i++)
        expect[i] = cos(t[i]);

    spida::SpidaRVT spi(grid);
    spi.dT(in,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);

    std::vector<dcmplx> out1(grid.getNst());
    spida::FFTRVT tr(grid);
    tr.T_To_ST(in,out1);
}