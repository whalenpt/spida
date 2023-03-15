
#include <gtest/gtest.h>
#include <fstream>
#include <complex>
#include <cmath>
#include <cassert>
#include <random>
#include <pwutils/pwmath.hpp>
#include "kiss_fft.h"
#include "kiss_fftr.h"


constexpr auto PI = 3.141592653589793238462643383279502884197;
using dcmplx = std::complex<double>;
constexpr dcmplx ii (0.0,1.0);

void createGrids(unsigned nx,double xmin,double xmax,std::vector<double>& x,std::vector<double>& kx)
{
    assert(x.size() == nx);
    assert(kx.size() == nx);
    double L = xmax - xmin;
    double dx = L/static_cast<double>(nx);
    for(unsigned i = 0; i < nx; i++) x[i] = xmin + i*dx; 
    double dsx = 2.0*PI/L;
    for(unsigned i = 0; i <= nx/2; i++)
        kx[i] = i*dsx;
    for(unsigned i = nx/2+1; i < nx; i++)
        kx[i] = -static_cast<double>(nx-i)*dsx;
}

/* Apply forward fft then inverse fft to a random vector and verify its the identity */
TEST(KISS_TEST,RANDOM)
{
    // Check that |in| = |ifft(fft(in))|, where ifft is the inverse fft
	unsigned N = 64;

    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<std::complex<double> > in(N);
    std::vector<std::complex<double> > out(N);
    for(unsigned int i = 0; i < N; i++)
        in[i] = distribution(generator);

    double in_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
    in_norm = sqrt(in_norm);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft_cfg cfg_reverse = kiss_fft_alloc(N,1,nullptr,nullptr);

    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft(cfg_reverse,reinterpret_cast<kiss_fft_cpx*>(out.data()),\
        reinterpret_cast<kiss_fft_cpx*>(in.data()));

    kiss_fft_free(cfg_forward);
    kiss_fft_free(cfg_reverse);

    for(unsigned i = 0; i < N; i++)
        in[i] = in[i]/static_cast<double>(N);

    double after_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
    after_norm = sqrt(after_norm);

    EXPECT_DOUBLE_EQ(in_norm,after_norm);
}

TEST(KISS_TEST,RANDOM_REAL)
{
    // Check that |in| = |irfft(rfft(in))| where irfft is inverse real-data fft 
	unsigned N = 64;

    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<double> in(N);
    std::vector<std::complex<double> > out(N/2+1);
    for(unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    auto in_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
    in_norm = sqrt(in_norm);

    kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
    kiss_fftr_cfg rcfg_reverse = kiss_fftr_alloc(N,1,nullptr,nullptr);

    kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fftri(rcfg_reverse,reinterpret_cast<kiss_fft_cpx*>(out.data()),\
            reinterpret_cast<kiss_fft_scalar*>(in.data()));

    kiss_fft_free(rcfg_forward);
    kiss_fft_free(rcfg_reverse);

    for(unsigned i = 0; i < N; i++)
        in[i] = in[i]/static_cast<double>(N);

    auto after_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const double& val) {return sum + pow(val,2);});
    after_norm = sqrt(after_norm);
    EXPECT_DOUBLE_EQ(in_norm,after_norm);
}

/* Inputs to FFT are all zeros, check that the sum of the output is also zero */
TEST(KISS_TEST,ZEROS)
{
	unsigned N = 64;
    std::vector<std::complex<double> > in(N, 0.);
    std::vector<std::complex<double> > out(N, 0.);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);

    auto zero_sum = std::accumulate(out.cbegin(),out.cend(),0.0,
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
	EXPECT_DOUBLE_EQ(zero_sum,0.0);
}

/* Inputs to FFT are all zeros, check that the sum of the output is also zero */
TEST(KISS_TEST,ZEROS_REAL)
{
	int N = 64;

    std::vector<double> in(N,0.0);
    std::vector<std::complex<double> >out(N/2+1);

    kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
    kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(rcfg_forward);

    auto zero_sum = std::accumulate(out.cbegin(),out.cend(),0.0,
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
	EXPECT_DOUBLE_EQ(zero_sum,0.0);
}

/* Inputs to FFT are all ones, check that the spectrum is centered at the zero DC component.*/ 
TEST(KISS_TEST,ONES)
{
	int N = 64;
    std::vector<std::complex<double> > in(N,1.0);
    std::vector<std::complex<double> > out(N);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);
	EXPECT_DOUBLE_EQ(out[0].real(),static_cast<double>(N));
}

/* Inputs to FFT are all ones, check that the spectrum is centered at the zero DC component.*/ 
TEST(KISS_TEST,ONES_REAL)
{
	int N = 64;
    std::vector<double> in(N,1.0);
    std::vector<std::complex<double> > out(N/2+1);

    kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
    kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(rcfg_forward);
	EXPECT_DOUBLE_EQ(out[0].real(),static_cast<double>(N));
}

/* Inputs to FFT are alternating 1,-1 */
TEST(KISS_TEST,ALTERNATING_ONE_NEGONE)
{
	unsigned N = 64;
    std::vector<std::complex<double> > in(N,1.0);
    std::vector<std::complex<double> > out(N);
	for(unsigned i = 1; i < in.size(); i+=2) in[i] = -1.0;

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));
	kiss_fft_free(cfg_forward);

	EXPECT_EQ(out[N/2].real(),static_cast<double>(N));
}

/* Inputs to FFT are alternating 1,-1 */
TEST(KISS_TEST,ALTERNATING_ONE_NEGONE_REAL)
{
	int N = 64;
    std::vector<double> in(N,1.0);
    std::vector<std::complex<double> > out(N/2+1);
	for(size_t i = 1; i < in.size(); i+=2) in[i] = -1.0;

	kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));
	kiss_fft_free(rcfg_forward);

	EXPECT_EQ(out.back().real(),static_cast<double>(N));
}

/* Inputs to FFT are imaginary e^(8*j*2*pi*i/N) for i = 0,1,2, ...,N-1 */
TEST(KISS_TEST,EXP_IMAG_TRIG_WAVE)
{
	unsigned N = 32;

    double xmin = 0.0;
    double xmax = 1.0;
    std::vector<double> x(N);
    std::vector<double> kx(N);
    createGrids(N,xmin,xmax,x,kx);
    std::vector<std::complex<double> > in(N, 0.);
    std::vector<std::complex<double> > out(N, 0.);
    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);

    for(unsigned j = 0; j < 8; j++){
        for(unsigned i = 0; i < N; i++)
            in[i] = exp(ii*2.0*PI*static_cast<double>(j)*x[i]);
        kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
                reinterpret_cast<kiss_fft_cpx*>(out.data()));
	    EXPECT_DOUBLE_EQ(out[j].real(),static_cast<double>(N));
    }

    kiss_fft_free(cfg_forward);
}


/* Inputs to FFT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 */
TEST(KISS_TEST,COS_TRIG_WAVE)
{
	unsigned N = 64;
	std::vector<std::complex<double>> in(N,0.0);
    for(unsigned i = 0; i < N; i++){
        in[i] = cos(8.0*2.0*PI*i/static_cast<double>(N));
    }
	std::vector<std::complex<double>> out(N);
	std::vector<std::complex<double>> out2(N);

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);

	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    for(unsigned i = 0; i < N; i++){
        in[i] = cos((43.0/7.0)*2.0*PI*i/static_cast<double>(N));
    }
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out2.data()));

    kiss_fft_free(cfg_forward);

    // Check cosine wave with frequency that lands on discretized grid
	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N)/2);

    // Check cosine wave with frequency that is in between grid frequencies
	// Compare results to python numpy rfft results for several values
	EXPECT_NEAR(out2[0].real(),1.445134428062532,1e-8);
	EXPECT_NEAR(out2[3].imag(),0.4098236995245791,1e-8);
	EXPECT_NEAR(out2[N/2].real(),0.06667213305344744,1e-8);
}

/* Inputs to FFT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 */
TEST(KISS_TEST,COS_TRIG_WAVE_REAL)
{
	unsigned N = 64;
	std::vector<double> in(N,0.0);
    for(unsigned int i = 0; i < N; i++){
        in[i] = cos(8.0*2.0*PI*i/static_cast<double>(N));
    }
	std::vector<std::complex<double>> out(N/2+1);

	kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    // Check cosine wave with frequency that lands on discretized grid
	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N)/2);

	std::vector<std::complex<double>> out2(N/2+1);
    for(unsigned i = 0; i < 64; i++){
        in[i] = cos((43.0/7.0)*2.0*i/static_cast<double>(N));
    }
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out2.data()));

    // Check cosine wave with frequency that is in between grid frequencies
	// Compare results to python numpy rfft results for several values
	EXPECT_NEAR(out2[0].real(),-1.4189089192588156,1e-8);
	EXPECT_NEAR(out2[3].imag(),-0.229996893956892,1e-8);
	EXPECT_NEAR(out2[32].real(),0.032896917953961546,1e-8);
	EXPECT_NEAR(out2[32].imag(),0.0,1e-8);
    kiss_fft_free(rcfg_forward);
}

TEST(KISS_TEST,GAUSS)
{
	unsigned N = 64;
    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> phase_adj(N,0.0);
    std::vector<dcmplx> out_phase_adj(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    double xmin = -6.0;
    double xmax = 6.0;
    double L = xmax-xmin;
    std::vector<double> x(N);
    std::vector<double> kx(N);
    createGrids(N,xmin,xmax,x,kx);

    double alpha = 2.0;
    for(size_t i = 0; i < x.size(); i++){
        in[i] = exp(-alpha*pow(x[i],2));
    }
    for(size_t i = 0; i < kx.size(); i++){
        expect[i] = (1.0/L)*sqrt(PI/alpha)*exp(-pow(kx[i],2)/(4.0*alpha));
    }

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    for(auto& item : out){
        item /= static_cast<double>(N);
    }

	for(size_t i = 0; i < out.size(); i++){
	    out_phase_adj[i] = out[i]*exp(ii*kx[i]*xmin);
    }

    EXPECT_LT(pw::relative_error(expect,out_phase_adj),1e-5);
}


TEST(KISS_TEST,SECH)
{
	unsigned N = 64;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    // make range large enough such that sech is near zero at xmin and xmax (FFT is periodic)
    double xmin = -6.0;
    double xmax = 6.0;
    double L = xmax-xmin;
    std::vector<double> x(N);
    std::vector<double> kx(N);
    createGrids(N,xmin,xmax,x,kx);
    double a = 2.0;
    // fill sech
    for(size_t i = 0; i < x.size(); i++){
        in[i] = 1.0/cosh(a*x[i]);
    }
    // spectrum of sech function
    for(size_t i = 0; i < kx.size(); i++){
        expect[i] = (1.0/L)*(PI/a)/cosh(PI*kx[i]/(2.0*a));
    }

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    // divide by FFT multiplier
    for(auto& item : out){
        item /= static_cast<double>(N);
    }

	for(size_t i = 0; i < out.size(); i++){
	    out[i] = out[i]*exp(ii*kx[i]*xmin);
    }

    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

// Test exp_imag_trig_wave 
TEST(KISS_TEST,EXP_IMAG_TRIG_WAVE2)
{
	unsigned N = 32;

    double xmin = 0.0;
    double xmax = 4.0;
    std::vector<double> x(N);
    std::vector<double> kx(N);
    createGrids(N,xmin,xmax,x,kx);

    std::vector<std::complex<double> > in(N, 0.);
    for(unsigned i = 0; i < N; i++)
        in[i] = exp(4.0*ii*2.0*PI*x[i]);
    std::vector<std::complex<double> > out(N, 0.);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);
	EXPECT_DOUBLE_EQ(out[16].real(),static_cast<double>(N));
}
