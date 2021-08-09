
#include <gtest/gtest.h>
#include <spida/constants.h>
#include <spida/grid/uniformT.h>
#include <spida/shape/shapeT.h>
#include <spida/grid/chebX.h>
#include <spida/transform/chebX.h>
#include <nayukidct/FastDctLee.hpp>
#include <nayukidct/FastDctFft.hpp>
#include <spida/spidaChebX.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <random>
#include <fftw3.h>


TEST(DCT_TEST,INVERSES)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);
    for(unsigned int i = 0; i < N; i++)
        in[i] = distribution(generator);

    spida::ChebTransformX tr(spida::ChebRootGridX{N});
    tr.X_To_SX(in,out);
    tr.SX_To_X(out,expect);

    dataio.writeFile("dct_fastlee_check.dat",in,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

TEST(DCT_TEST,TRANSFORM_TEST)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    spida::ChebRootGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < N; i++)
        in[i] = exp(x[i]);

    spida::ChebTransformX tr(grid);
    std::vector<double> out(N);
    tr.X_To_SX(in,out);

    std::vector<double> out2(N);
    std::reverse_copy(std::cbegin(in),std::cend(in),std::begin(out2));
    FastDctLee::transform(out2);
    for(auto& item : out2)
        item /= (N/2.0);

    dataio.writeFile("dct_transform_test.dat",out,out2);
    EXPECT_LT(pw::relative_error(out,out2),1e-6);
}

TEST(DCT_TEST,DERIVATIVE_EXP)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    std::vector<double> out(N,0.0);

    spida::ChebRootGridX grid(N,-1,5);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(x[i]);

    spida::ChebX spidaX(grid);
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(in,out),1e-6);
}

TEST(DCT_TEST,COMP_SPECTRA)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    std::vector<double> out(N);

    spida::ChebRootGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(x[i]);
    std::copy(std::cbegin(in),std::cend(in),std::begin(out));
    FastDctLee::transform(out);

    std::vector<double> out2(N);
    fftw_plan forward_plan = fftw_plan_r2r_1d(N,in.data(),out2.data(),FFTW_REDFT10,FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);

    for(auto& item : out)
        item *= 2.0;

    EXPECT_LT(pw::relative_error(out,out2),1e-6);
    dataio.writeFile("spectra_compare.dat",out,out2);
}

TEST(DCT_TEST,NAYUKI_ROOT_DERS)
{
	int N = 16;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    spida::ChebRootGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(x[i]);

    std::vector<double> alpha(N);
    std::copy(std::cbegin(in),std::cend(in),std::begin(alpha));
    FastDctLee::transform(alpha);

    dataio.writeFile("nayuki_alphaval.dat",alpha);
    for(auto& item : alpha)
        item /= static_cast<double>(-N/2.0);

    std::vector<double> beta(N);
    beta[N-1] = 0.0;
    beta[N-2] =  2.0*(N-1)*alpha[N-1]; 
    for(int k = N-2; k > 0; k--)
        beta[k-1] = beta[k+1] + 2.0*k*alpha[k];

    FastDctLee::inverseTransform(beta);
    dataio.writeFile("nayuki_der_compare.dat",in,beta);
    EXPECT_LT(pw::relative_error(in,beta),1e-6);
}

TEST(DCT_TEST,DCT_EXP_DER_TEST)
{
	int N = 16;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    spida::ChebRootGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(x[i]);

    std::vector<double> alpha(N);
    std::copy(std::cbegin(in),std::cend(in),std::begin(alpha));
    FastDctLee::transform(alpha);

    dataio.writeFile("nayuki_alphaval.dat",alpha);
    for(auto& item : alpha)
        item /= static_cast<double>(-N/2.0);

    std::vector<double> beta(N);
    beta[N-1] = 0.0;
    beta[N-2] =  2.0*(N-1)*alpha[N-1]; 
    for(int k = N-2; k > 0; k--)
        beta[k-1] = beta[k+1] + 2.0*k*alpha[k];

    FastDctLee::inverseTransform(beta);
    dataio.writeFile("nayuki_der_compare.dat",in,beta);
    EXPECT_LT(pw::relative_error(in,beta),1e-6);
}





