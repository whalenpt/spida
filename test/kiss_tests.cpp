
#include <gtest/gtest.h>
#include <fstream>
#include <complex>
#include <cmath>
#include <random>
#include <pwutils/report/dataio.hpp>
#include <spida/helper/constants.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"

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
	int N = 64;
    pw::DataIO dataio("outfolder");
    std::vector<std::complex<double> > in(N, 0.);
    for(unsigned int i = 0; i < N; i++)
        in[i] = exp(8.0*spida::ii*2.0*spida::PI*static_cast<double>(i)/static_cast<double>(N));
    std::vector<std::complex<double> > out(N, 0.);

    kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);
    kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
            reinterpret_cast<kiss_fft_cpx*>(out.data()));
    kiss_fft_free(cfg_forward);
    dataio.writeFile("kiss_exp_trig.dat",in,out);
	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N));
}

/* Inputs to FFT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 */
TEST(KISS_TEST,COS_TRIG_WAVE)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

	std::vector<std::complex<double>> in(N,0.0);
    for(auto i = 0; i < N; i++)
        in[i] = cos(8.0*2.0*spida::PI*i/static_cast<double>(N));
	std::vector<std::complex<double>> out(N);
	std::vector<std::complex<double>> out2(N);

	kiss_fft_cfg cfg_forward = kiss_fft_alloc(N,0,nullptr,nullptr);

	kiss_fft(cfg_forward,reinterpret_cast<kiss_fft_cpx*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));

    for(unsigned int i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*2.0*spida::PI*i/static_cast<double>(N));
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
        in[i] = cos(8.0*2.0*spida::PI*i/static_cast<double>(N));
	std::vector<std::complex<double>> out(N/2+1);
	std::vector<std::complex<double>> out2(N/2+1);

	kiss_fftr_cfg rcfg_forward = kiss_fftr_alloc(N,0,nullptr,nullptr);
	kiss_fftr(rcfg_forward,reinterpret_cast<kiss_fft_scalar*>(in.data()),\
	        reinterpret_cast<kiss_fft_cpx*>(out.data()));
    for(auto i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*2.0*spida::PI*i/static_cast<double>(N));
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










