/*
MIT License

Copyright (c) 2021 Patrick Whalen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <gtest/gtest.h>
#include <nayukidct/FastDctLee.hpp>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>

// Check that transform followed by inverse transform is identity
TEST(NAYUKI_TEST,RANDOM)
{
	int N = 64;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<double> in(N);
    std::vector<double> out(N);
    for(auto i = 0; i < N; i++)
        in[i] = distribution(generator);

    std::copy(in.cbegin(),in.cend(),out.begin());
    FastDctLee::transform(out);
    FastDctLee::inverseTransform(out);

    for(auto& item : out)
        item *= 2.0/static_cast<double>(N);

    double diff = 0.0;
    for(auto i = 0; i < in.size(); i++)
        diff += abs(in[i] - out[i]);
    EXPECT_NEAR(diff,0.0,1e-12);
}

// Inputs to DCT are all zeros, check that the sum of the output is also zero
TEST(NAYUKI_TEST,ZEROS)
{
	int N = 64;

    std::vector<double> in(N, 0.);
    std::vector<double> out(N);
    std::copy(in.cbegin(),in.cend(),out.begin());
    FastDctLee::transform(out);
	double sum = 0.0;
	for(auto& val : out){
		sum += abs(val); 
	}
	EXPECT_EQ(sum,0.0);
}

// Inputs to DCT are all ones, check that the spectrum is centered at the zero DC component. 
TEST(NAYUKI_TEST,ONES)
{
	int N = 64;
    std::vector<double> in(N, 1.);
    std::vector<double> out(N);
    std::copy(in.cbegin(),in.cend(),out.begin());
    FastDctLee::transform(out);
	EXPECT_EQ(out[0],static_cast<double>(N));
}

// Inputs to DCT are alternating 1,-1 
TEST(NAYUKI_TEST,ALTERNATING_ONE_NEGONE)
{
	int N = 64;
    std::vector<double> in(N,1.0);
	for(auto i = 1; i < in.size(); i+=2)
		in[i] = -1.0;
    std::vector<double> out(N);
    std::copy(in.cbegin(),in.cend(),out.begin());
    FastDctLee::transform(out);

    for(auto i = 0; i < N; i++)
        out[i] = 2*out[i]/N;

    // Results compared to fftw_plan_r2r  (FFTW)
	EXPECT_NEAR(out[0],0.0,1e-6);
	EXPECT_NEAR(out[1],3.125941475e-02,1e-6);
	EXPECT_NEAR(out[N-1],1.27336738,1e-6);
}

// Inputs to DCT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 
TEST(NAYUKI_TEST,COS_TRIG_WAVE)
{
    double pi = 3.141592653589793238462643383279502884;
	int N = 64;
    std::vector<double> in(N);
    for(auto i = 0; i < N; i++)
        in[i] = cos(8*pi*(i+0.5)/static_cast<double>(N));
    std::vector<double> out(N);
    std::copy(in.cbegin(),in.cend(),out.begin());
    FastDctLee::transform(out);
	// Trig wave consists of sines and cosines at frequency bin 8
	EXPECT_DOUBLE_EQ(out[8],static_cast<double>(N)/2);
}

// Inputs to DCT are a cosine wave
TEST(NAYUKI_TEST,COS_TRIG_WAVE2)
{
    double pi = 3.141592653589793238462643383279502884;
	int N = 64;
    std::vector<double> in(N);
    for(auto i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*pi*(i+0.5)/static_cast<double>(N));
    std::vector<double> out(N);
    std::copy(in.cbegin(),in.cend(),out.begin());
    FastDctLee::transform(out);
    for(auto i = 0; i < N; i++)
        out[i] = 2*out[i]/static_cast<double>(N);
    // Results compared to fftw_plan_r2r  (FFTW)
	EXPECT_NEAR(out[5],-1.3341195346634438e-01,1e-6);
	EXPECT_NEAR(out[6],9.7831222281714325e-01,1e-6);
	EXPECT_NEAR(out[7],1.5044681031652760e-01,1e-6);
}

// Check that transform followed by inverse transform is identity
TEST(PLAN_TEST,RANDOM)
{
	int N = 64;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);
    for(auto i = 0; i < N; i++)
        in[i] = distribution(generator);

    FastDctPlan plan_forward(in,out,dct_kind::DCTII);
    plan_forward.execute();
    FastDctPlan plan_backward(out,expect,dct_kind::DCTIII);
    plan_backward.execute();

    for(auto& item : expect)
        item *= 2.0/static_cast<double>(N);

    double diff = 0.0;
    for(auto i = 0; i < in.size(); i++)
        diff += abs(in[i] - expect[i]);
    EXPECT_NEAR(diff,0.0,1e-12);
}


// Inputs to DCT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 
TEST(PLAN_TEST,COS_TRIG_WAVE)
{
    double pi = 3.141592653589793238462643383279502884;
	int N = 64;
    std::vector<double> in(N);
    for(auto i = 0; i < N; i++)
        in[i] = cos(8*pi*(i+0.5)/static_cast<double>(N));
    std::vector<double> out(N);
    FastDctPlan plan_forward(in,out,dct_kind::DCTII);
    plan_forward.execute();
	// Trig wave consists of sines and cosines at frequency bin 8
	EXPECT_DOUBLE_EQ(out[8],static_cast<double>(N)/2);
}

TEST(NAYUKI_TEST,EXP_DER_TEST)
{
	int N = 64;
    double pi = 3.141592653589793238462643383279502884;

    // Set up Chebyshev Root Grid on [a,b]
    std::vector<double> x(N);
    double a = -2.0;
    double b = 4.0;
    double dc = pi/static_cast<double>(N);
    double L = b - a;
    //for(auto j = 0; j < N; j++) x[N-j-1] = (cos((j+0.5)*dc)*L+b+a)/2.0; 
    for(auto j = 0; j < N; j++) x[j] = (cos((j+0.5)*dc)*L+b+a)/2.0; 

    // Setup exponential input
    std::vector<double> in(N);
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(2*x[i]);

    std::vector<double> alpha(N);
    std::copy(in.cbegin(),in.cend(),alpha.begin());
    FastDctLee::transform(alpha);

    for(auto& item : alpha)
        item /= static_cast<double>(N/2.0);

    // Take spectral derivative
    std::vector<double> beta(N);
    double sf = 2.0/L;
    beta[N-1] = 0.0;
    beta[N-2] =  2.0*(N-1)*sf*alpha[N-1]; 
    for(int k = N-2; k > 0; k--)
        beta[k-1] = beta[k+1] + 2.0*sf*k*alpha[k];

    FastDctLee::inverseTransform(beta);

    double diff = 0.0;
    for(auto i = 0; i < in.size(); i++)
        diff += abs(2.0*in[i] - beta[i]);
    diff /= N;
    EXPECT_NEAR(diff,0.0,1e-10);
}


#include <cassert>
#include <cmath>
TEST(NAYUKI_TEST,EXP_DER_TEST2)
{
	int N = 64;
    double pi = 3.141592653589793238462643383279502884;

    // Set up Chebyshev Root Grid on [a,b]
    std::vector<double> x(N);
    double dc = pi/static_cast<double>(N);
    for(auto j = 0; j < N; j++) x[j] = cos((j+0.5)*dc); 

    // Setup exponential input
    std::vector<double> in(N);
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(4.0*x[i]);

    std::vector<double> alpha(N);
    std::copy(in.cbegin(),in.cend(),alpha.begin());
    FastDctLee::transform(alpha);
    for(auto& item : alpha)
        item /= static_cast<double>(N/2.0);

    // Take spectral derivative
    std::vector<double> beta(N);
    beta[N-1] = 0.0;
    beta[N-2] =  2.0*(N-1)*alpha[N-1]; 
    for(int k = N-2; k > 0; k--)
        beta[k-1] = beta[k+1] + 2.0*k*alpha[k];

    FastDctLee::inverseTransform(beta);

    double diff = 0.0;
    for(auto i = 0; i < in.size(); i++)
        diff += abs(4.0*in[i] - beta[i]);
    diff /= N;
    EXPECT_NEAR(diff,0.0,1e-10);
}





