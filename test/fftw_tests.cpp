
#include <gtest/gtest.h>
#include <memory>
#include <fstream>
#include <complex>
#include <cmath>
#include <random>
#include <pwutils/report/dataio.hpp>
#include <spida/propagator/propagator.h>
#include <spida/constants.h>
#include "fftw3.h"

// Apply forward fft then inverse fft to a random vector and verify its the identity 

TEST(FFTW_TEST,RANDOM)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<std::complex<double> > in(N);
    std::vector<std::complex<double> > out(N, 0.);
    for(auto i = 0; i < N; i++)
        in[i] = distribution(generator);

    double in_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
    in_norm = sqrt(in_norm);


    fftw_plan plan_forward = fftw_plan_dft_1d(N,reinterpret_cast<fftw_complex*>(in.data()),\
            reinterpret_cast<fftw_complex*>(out.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_reverse = fftw_plan_dft_1d(N,reinterpret_cast<fftw_complex*>(out.data()),\
            reinterpret_cast<fftw_complex*>(in.data()), FFTW_BACKWARD, FFTW_ESTIMATE);

    dataio.writeFile("random_before.dat",in);
    fftw_execute(plan_forward);
    fftw_execute(plan_reverse);

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_reverse);

    for(auto i = 0; i < N; i++)
        in[i] = in[i]/static_cast<double>(N);

    dataio.writeFile("random_after.dat",in);

    double after_norm = std::accumulate(in.cbegin(),in.cend(),0.0,\
            [](double const& sum,const std::complex<double>& val) {return sum + std::norm(val);});
    after_norm = sqrt(after_norm);
    dataio.appendFile("random_norms.dat",in_norm,after_norm);

    // Check that |in| = |ifft(fft(in))|, where ifft is the inverse fft
    EXPECT_DOUBLE_EQ(in_norm,after_norm);


    std::vector<double> in_rfft(N);
    std::vector<std::complex<double> > out_rfft(N/2+1);
    for(auto i = 0; i < N; i++)
        in_rfft[i] = in[i].real();


    fftw_plan rplan_forward = fftw_plan_dft_r2c_1d(N,in_rfft.data(),\
            reinterpret_cast<fftw_complex*>(out_rfft.data()), FFTW_ESTIMATE);
    fftw_plan rplan_reverse = fftw_plan_dft_c2r_1d(N,\
            reinterpret_cast<fftw_complex*>(out_rfft.data()),in_rfft.data(), FFTW_ESTIMATE);

    fftw_execute(rplan_forward);
    fftw_execute(rplan_reverse);
    fftw_destroy_plan(rplan_forward);
    fftw_destroy_plan(rplan_reverse);

    for(auto i = 0; i < N; i++)
        in_rfft[i] = in_rfft[i]/static_cast<double>(N);

    double after_norm2 = std::accumulate(in_rfft.cbegin(),in_rfft.cend(),0.0,\
            [](double const& sum,const double& val) {return sum + pow(val,2);});
    after_norm2 = sqrt(after_norm2);
    dataio.appendFile("random_norms.dat",in_norm,after_norm2);

    // Check that |in| = |irfft(rfft(in))| where irfft is inverse real-data fft 
    EXPECT_DOUBLE_EQ(in_norm,after_norm2);
}



// Inputs to FFT are all zeros, check that the sum of the output is also zero
TEST(FFTW_TEST,ZEROS)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<std::complex<double> > in(N, 0.);
    std::vector<std::complex<double> > out(N, 0.);
    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&in[0]),\
            reinterpret_cast<fftw_complex*>(&out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    dataio.writeFile("zeros_fft.dat",in,out);
	double sum = 0.0;
	for(auto& val : out){
		sum += abs(val); 
	}
	EXPECT_EQ(sum,0.0);
}

// Inputs to FFT are all ones, check that the spectrum is centered at the zero DC component. 
TEST(FFTW_TEST,ONES)
{
	int N = 64;
    pw::DataIO dataio("outfolder");
    std::vector<std::complex<double> > in(N, 1.);
    std::vector<std::complex<double> > out(N, 0.);

    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&in[0]),\
            reinterpret_cast<fftw_complex*>(&out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    dataio.writeFile("ones_fft.dat",in,out);
	EXPECT_EQ(out[0].real(),static_cast<double>(N));
}

// Inputs to FFT are alternating 1,-1 
TEST(FFTW_TEST,ALTERNATING_ONE_NEGONE)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N,1.0);
	for(auto i = 1; i < in.size(); i+=2)
		in[i] = -1.0;

	std::vector<std::complex<double>> out(N/2+1, 0.);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N,&in[0],reinterpret_cast<fftw_complex*>(&out[0]), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    dataio.writeFile("alternating_one_neg_one_fft.dat",out);
	EXPECT_EQ(out.back().real(),static_cast<double>(N));
}

// Inputs to FFT are imaginary e^(8*j*2*pi*i/N) for i = 0,1,2, ...,N-1 
TEST(FFTW_TEST,EXP_IMAG_TRIG_WAVE)
{
	int N = 64;
    pw::DataIO dataio("outfolder");
    std::vector<std::complex<double> > in(N, 0.);
    for(auto i = 0; i < N; i++)
        in[i] = exp(8.0*spida::ii*2.0*spida::PI*static_cast<double>(i)/static_cast<double>(N));
    std::vector<std::complex<double> > out(N, 0.);

    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&in[0]),\
            reinterpret_cast<fftw_complex*>(&out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    dataio.writeFile("exp_imag_trig_fft.dat",in,out);
	// Trig wave consists of sines and cosines at frequency bin 8
	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N));
}

// Inputs to FFT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 
TEST(FFTW_TEST,COS_TRIG_WAVE)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N,0.0);
    for(auto i = 0; i < N; i++)
        in[i] = cos(8.0*2.0*spida::PI*i/static_cast<double>(N));

	std::vector<std::complex<double>> out(N/2+1, 0.);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N,in.data(),reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    dataio.writeFile("cos_fft.dat",out);
	// Trig wave consists of sines and cosines at frequency bin 8
	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N)/2);

    for(auto i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*2.0*spida::PI*i/static_cast<double>(N));
    fftw_plan plan2 = fftw_plan_dft_r2c_1d(N,in.data(),reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE);
    fftw_execute(plan2);
    fftw_destroy_plan(plan2);

    dataio.writeFile("cos_fft2.dat",out);
	// Compare results to python numpy rfft results for several values
	EXPECT_NEAR(out[0].real(),1.44513443,1e-6);
	EXPECT_NEAR(out[3].imag(),0.4098237,1e-6);
	EXPECT_NEAR(out[N/2].real(),6.66721331e-2,1e-6);
}










