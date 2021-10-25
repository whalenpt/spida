
#include <gtest/gtest.h>
#include <fstream>
#include <complex>
#include <cmath>
#include <cassert>
#include <random>
#include <pwutils/report/dataio.hpp>
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
    for(auto i = 0; i < nx; i++) x[i] = xmin + i*dx; 
    double dsx = 2.0*PI/L;
    for(auto i = 0; i <= nx/2; i++)
        kx[i] = i*dsx;
    for(auto i = nx/2+1; i < nx; i++)
        kx[i] = -static_cast<double>(nx-i)*dsx;
}

/* Apply forward fft then inverse fft to a random vector and verify its the identity */
TEST(KISS_TEST,RANDOM)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<std::complex<double> > in(N);
    std::vector<std::complex<double> > out(N);
    for(unsigned int i = 0; i < N; i++)
        in[i] = distribution(generator);
    dataio.writeFile("kiss_random_before.dat",in);

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

    for(auto i = 0; i < N; i++)
        in[i] = in[i]/static_cast<double>(N);

    dataio.writeFile("kiss_random_after.dat",in);

    double after_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
    after_norm = sqrt(after_norm);
    dataio.appendFile("kiss_random_norms.dat",in_norm,after_norm);

    // Check that |in| = |ifft(fft(in))|, where ifft is the inverse fft
    EXPECT_DOUBLE_EQ(in_norm,after_norm);
}

TEST(KISS_TEST,RANDOM_REAL)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<double> in(N);
    std::vector<std::complex<double> > out(N/2+1);
    for(unsigned int i = 0; i < N; i++)
        in[i] = distribution(generator);
    dataio.writeFile("kiss_real_random_before.dat",in);

    double in_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
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

    for(auto i = 0; i < N; i++)
        in[i] = in[i]/static_cast<double>(N);

    dataio.writeFile("kiss_real_random_after.dat",in);

    double after_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const double& val) {return sum + pow(val,2);});
    after_norm = sqrt(after_norm);
    dataio.appendFile("kiss_real_random_norms.dat",in_norm,after_norm);

    // Check that |in| = |irfft(rfft(in))| where irfft is inverse real-data fft 
    EXPECT_DOUBLE_EQ(in_norm,after_norm);
}

/* Inputs to FFT are all zeros, check that the sum of the output is also zero */
TEST(KISS_TEST,ZEROS)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<std::complex<double> > in(N, 0.);
    std::vector<std::complex<double> > out(N, 0.);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);

    dataio.writeFile("kiss_zeros.dat",in,out);
    double sum = std::accumulate(out.cbegin(),out.cend(),0.0,
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
	EXPECT_DOUBLE_EQ(sum,0.0);
}

/* Inputs to FFT are all zeros, check that the sum of the output is also zero */
TEST(KISS_TEST,ZEROS_REAL)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N,0.0);
    std::vector<std::complex<double> >out(N/2+1);

    kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
    kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(rcfg_forward);

    double sum = std::accumulate(out.cbegin(),out.cend(),0.0,
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
	EXPECT_DOUBLE_EQ(sum,0.0);
}

/* Inputs to FFT are all ones, check that the spectrum is centered at the zero DC component.*/ 
TEST(KISS_TEST,ONES)
{
	int N = 64;
    pw::DataIO dataio("outfolder");
    std::vector<std::complex<double> > in(N,1.0);
    std::vector<std::complex<double> > out(N);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);
    dataio.writeFile("kiss_ones.dat",in,out);
	EXPECT_DOUBLE_EQ(out[0].real(),static_cast<double>(N));
}

/* Inputs to FFT are all ones, check that the spectrum is centered at the zero DC component.*/ 
TEST(KISS_TEST,ONES_REAL)
{
	int N = 64;
    pw::DataIO dataio("outfolder");
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
	int N = 64;
    pw::DataIO dataio("outfolder");
    std::vector<std::complex<double> > in(N,1.0);
    std::vector<std::complex<double> > out(N);
	for(auto i = 1; i < in.size(); i+=2) in[i] = -1.0;

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));
	kiss_fft_free(cfg_forward);

    dataio.writeFile("kiss_alternating_one_neg_one.dat",in,out);
	EXPECT_EQ(out[N/2].real(),static_cast<double>(N));
}

/* Inputs to FFT are alternating 1,-1 */
TEST(KISS_TEST,ALTERNATING_ONE_NEGONE_REAL)
{
	int N = 64;
    pw::DataIO dataio("outfolder");
    std::vector<double> in(N,1.0);
    std::vector<std::complex<double> > out(N/2+1);
	for(auto i = 1; i < in.size(); i+=2) in[i] = -1.0;

	kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));
	kiss_fft_free(rcfg_forward);

	EXPECT_EQ(out.back().real(),static_cast<double>(N));
}

/* Inputs to FFT are imaginary e^(8*j*2*pi*i/N) for i = 0,1,2, ...,N-1 */
TEST(KISS_TEST,EXP_IMAG_TRIG_WAVE)
{
	int N = 32;

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
	int N = 64;
    pw::DataIO dataio("outfolder");

	std::vector<std::complex<double>> in(N,0.0);
    for(auto i = 0; i < N; i++)
        in[i] = cos(8.0*2.0*PI*i/static_cast<double>(N));
	std::vector<std::complex<double>> out(N);
	std::vector<std::complex<double>> out2(N);

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);

	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    for(unsigned int i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*2.0*PI*i/static_cast<double>(N));
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out2.data()));

    kiss_fft_free(cfg_forward);

    dataio.writeFile("kiss_cos_fft.dat",out);
    dataio.writeFile("kiss_cos_fft2.dat",out2);

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
	int N = 64;
    pw::DataIO dataio("outfolder");

	std::vector<double> in(N,0.0);
    for(unsigned int i = 0; i < N; i++)
        in[i] = cos(8.0*2.0*PI*i/static_cast<double>(N));
	std::vector<std::complex<double>> out(N/2+1);
	std::vector<std::complex<double>> out2(N/2+1);

	kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    for(auto i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*2.0*i/static_cast<double>(N));
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out2.data()));

    kiss_fft_free(rcfg_forward);

    // Check cosine wave with frequency that lands on discretized grid
	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N)/2);
    // Check cosine wave with frequency that is in between grid frequencies
	// Compare results to python numpy rfft results for several values
	EXPECT_NEAR(out2[0].real(),1.445134428062532,1e-8);
	EXPECT_NEAR(out2[3].imag(),0.4098236995245791,1e-8);
	EXPECT_NEAR(out2[N/2].real(),0.06667213305344744,1e-8);
}

TEST(KISS_TEST,GAUSS)
{
	unsigned N = 64;
    pw::DataIO dataio("outfolder");
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
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(-alpha*pow(x[i],2));
    for(auto i = 0; i < kx.size(); i++)
        expect[i] = (1.0/L)*sqrt(PI/alpha)*exp(-pow(kx[i],2)/(4.0*alpha));

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    for(auto& item : out)
        item /= static_cast<double>(N);

	for(auto i = 0; i < out.size(); i++)
	    out_phase_adj[i] = out[i]*exp(ii*kx[i]*xmin);
	    //out_phase_adj[i] = out[i]*exp(2.0*PI*ii*kx[i]*xmin);

    dataio.writeFile("kissfft_gauss.dat",expect,out_phase_adj);
    EXPECT_LT(pw::relative_error(expect,out_phase_adj),1e-5);
}


TEST(KISS_TEST,SECH)
{
	unsigned N = 64;
    pw::DataIO dataio("outfolder");

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
    for(auto i = 0; i < x.size(); i++)
        in[i] = 1.0/cosh(a*x[i]);
    // spectrum of sech function
    for(auto i = 0; i < kx.size(); i++)
        expect[i] = (1.0/L)*(PI/a)/cosh(PI*kx[i]/(2.0*a));

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    // divide by FFT multiplier
    for(auto& item : out)
        item /= static_cast<double>(N);

	for(auto i = 0; i < out.size(); i++)
	    out[i] = out[i]*exp(ii*kx[i]*xmin);

    dataio.writeFile("kissfft_sech.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

// Test exp_imag_trig_wave 
TEST(KISS_TEST,EXP_IMAG_TRIG_WAVE2)
{
	int N = 32;
    pw::DataIO dataio("outfolder");

    double xmin = 0.0;
    double xmax = 4.0;
    std::vector<double> x(N);
    std::vector<double> kx(N);
    createGrids(N,xmin,xmax,x,kx);

    std::vector<std::complex<double> > in(N, 0.);
    for(unsigned int i = 0; i < N; i++)
        in[i] = exp(4.0*ii*2.0*PI*x[i]);
    std::vector<std::complex<double> > out(N, 0.);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);
    dataio.writeFile("kiss_exp_trig2.dat",kx,out);
	EXPECT_DOUBLE_EQ(out[16].real(),static_cast<double>(N));
}













