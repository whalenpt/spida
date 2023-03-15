
#include <gtest/gtest.h>
#include <pwutils/pwmath.hpp>
#include <nayukidct/FastDctLee.hpp>
#include <random>
#include <numeric>

constexpr auto PI = 3.141592653589793238462643383279502884197;

TEST(NAYUKI_TEST,RANDOM)
{
    // Check that |in| = |ifft(fft(in))|, where ifft is the inverse fft
	unsigned N = 64;
    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<double> in(N);
    std::vector<double> out(N);
    for(unsigned i = 0; i < N; i++)
        in[i] = distribution(generator);

    double in_norm = sqrt(std::accumulate(in.cbegin(),in.cend(),0));
    std::copy(std::cbegin(in),std::cend(in),std::begin(out));
    FastDctLee::transform(out);
    FastDctLee::inverseTransform(out);

    for(unsigned i = 0; i < N; i++)
        out[i] = 2*out[i]/static_cast<double>(N);

    auto after_norm = sqrt(std::accumulate(out.cbegin(),out.cend(),0));
    EXPECT_DOUBLE_EQ(in_norm,after_norm);
}

// Inputs to DCT are all zeros, check that the sum of the output is also zero
TEST(NAYUKI_TEST,ZEROS)
{
	unsigned N = 64;
    std::vector<double> in(N, 0.);
    std::vector<double> out(N, 0.);
    std::copy(std::begin(in),std::end(in),std::begin(out));
    FastDctLee::transform(out);
	double sum = 0.0;
	for(const auto& val : out){
		sum += abs(val); 
	}
	EXPECT_EQ(sum,0.0);
}

// Inputs to DCT are all ones, check that the spectrum is centered at the zero DC component. 
TEST(NAYUKI_TEST,ONES)
{
	unsigned N = 64;
    std::vector<double> in(N, 1.);
    std::vector<double> out(N, 0.);
    std::copy(std::begin(in),std::end(in),std::begin(out));
    FastDctLee::transform(out);
	EXPECT_EQ(out[0],static_cast<double>(N));
}

// Inputs to DCT are alternating 1,-1 
TEST(NAYUKI_TEST,ALTERNATING_ONE_NEGONE)
{
	unsigned N = 64;
    std::vector<double> in(N,1.0);
	for(unsigned i = 1; i < in.size(); i+=2)
		in[i] = -1.0;
    std::vector<double> out(N);
    std::copy(std::cbegin(in),std::cend(in),std::begin(out));
    FastDctLee::transform(out);

    for(unsigned i = 0; i < N; i++){
        out[i] = 2*out[i]/N;
    }

	EXPECT_NEAR(out[0],0.0,1e-6);
	EXPECT_NEAR(out[1],3.125941475e-02,1e-6);
	EXPECT_NEAR(out[N-1],1.27336738,1e-6);
}

// Inputs to DCT are cos(8*2*pi*i/N) for i = 0,1,2, ...,N-1 
TEST(NAYUKI_TEST,COS_TRIG_WAVE)
{
	unsigned N = 64;
    std::vector<double> in(N,0.0);
    for(unsigned i = 0; i < N; i++)
        in[i] = cos(8*PI*(i+0.5)/static_cast<double>(N));
    std::vector<double> out(N);
    std::copy(std::cbegin(in),std::cend(in),std::begin(out));

    FastDctLee::transform(out);
	// Trig wave consists of sines and cosines at frequency bin 8
	EXPECT_DOUBLE_EQ(out[8],static_cast<double>(N)/2);


    for(unsigned i = 0; i < N; i++)
        in[i] = cos((43.0/7.0)*PI*(i+0.5)/static_cast<double>(N));
    std::copy(std::cbegin(in),std::cend(in),std::begin(out));
    FastDctLee::transform(out);
    for(unsigned i = 0; i < N; i++){
        out[i] = 2*out[i]/static_cast<double>(N);
    }

	// Compared results to fftw output
	EXPECT_NEAR(out[5],-1.3341195346634438e-01,1e-6);
	EXPECT_NEAR(out[6],9.7831222281714325e-01,1e-6);
	EXPECT_NEAR(out[7],1.5044681031652760e-01,1e-6);
}