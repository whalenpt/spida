
#include <gtest/gtest.h>
#include <spida/constants.h>
#include <spida/grid/uniformT.h>
#include <spida/shape/shapeT.h>
#include <spida/grid/chebX.h>
#include <spida/transform/chebFFTWX.h>
#include <spida/spidaFFTWX.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <random>
#include <fftw3.h>


TEST(FDCT_TEST,INVERSES)
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
    dataio.writeFile("dct_random_before.dat",in);

    spida::ChebExtremaGridX grid(N);
    spida::ChebTransformFFTWX trfft(grid);
    trfft.X_To_SX(in,out);
    trfft.SX_To_X(out,expect);
    dataio.writeFile("dct_fftw_check.dat",in,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}


TEST(FDCT_TEST,DERIVATIVE_EXP)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    std::vector<double> out(N,0.0);
    spida::ChebExtremaGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(x[i]);

    spida::ChebFFTWX spidaX(grid);
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(in,out),1e-6);
    dataio.writeFile("fdct_dexp.dat",in,out);

    std::vector<double> exact(N);
    for(auto i = 0; i < x.size(); i++){
        in[i] = exp(2*x[i]);
        exact[i] = 2*exp(2*x[i]);
    }
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);

    spida::ChebExtremaGridX grid2 = spida::ChebExtremaGridX(N,0,5);
    const std::vector<double> x2 = grid2.getX();
    for(auto i = 0; i < x2.size(); i++){
        in[i] = exp(-x2[i]);
        exact[i] = -exp(-x2[i]);
    }
    spida::ChebFFTWX spida2X(grid2);
    spida2X.dX(in,out);

    EXPECT_LT(pw::relative_error(out,exact),1e-6);
    dataio.writeFile("fdct_dexp2.dat",out,exact);
    
    spida::ChebExtremaGridX grid3 = spida::ChebExtremaGridX(N,-5,1);
    const std::vector<double> x3 = grid3.getX();
    for(auto i = 0; i < x3.size(); i++)
        in[i] = exp(x3[i]);
    spida::ChebFFTWX spida3X(grid3);
    spida3X.dX(in,out);
    dataio.writeFile("fdct_dexp3.dat",in,out);
    EXPECT_LT(pw::relative_error(in,out),1e-6);
}

TEST(FDCT_TEST,DERIVATIVE_COS)
{
	int N = 64;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    std::vector<double> exact(N);
    std::vector<double> out(N);
    spida::ChebExtremaGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++){
        in[i] = cos(x[i]);
        exact[i] = -sin(x[i]);
    }
    spida::ChebFFTWX spidaX(grid);
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);
    //dataio.writeFile("dct_dexp_extrema.dat",in,out);

    for(auto i = 0; i < x.size(); i++){
        in[i] = cos(3.6*x[i]);
        exact[i] = -3.6*sin(3.6*x[i]);
    }
    spidaX.dX(in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);


    spida::ChebExtremaGridX grid2(N,0,5);
    const std::vector<double> x2 = grid2.getX();
    for(auto i = 0; i < x2.size(); i++){
        in[i] = cos(3.6*x2[i]);
        exact[i] = -3.6*sin(3.6*x2[i]);
    }
    spida::ChebFFTWX spida2X(grid2);
    spida2X.dX(in,out);

    EXPECT_LT(pw::relative_error(out,exact),1e-6);
    dataio.writeFile("dct_dcos2_extrema.dat",out,exact);
    
    spida::ChebExtremaGridX grid3(N,-5,1);
    const std::vector<double> x3 = grid3.getX();
    for(auto i = 0; i < x3.size(); i++){
        in[i] = cos(3.6*x3[i]);
        exact[i] = -3.6*sin(3.6*x3[i]);
    }
    spida::ChebFFTWX spida3X(grid3);
    spida3X.dX(in,out);
    dataio.writeFile("dct_dcos3_extrema.dat",in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);
}

// Use ChebRootGrid with FFTW_REDFT10 (DCT-type II)  and FFTW_REDFT01 (DCT-type III)
TEST(FDCT_TEST,CHEBROOT_DERS)
{
	int N = 20;
    pw::DataIO dataio("outfolder");

    std::vector<double> in(N);
    spida::ChebRootGridX grid(N,-1,1);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(x[i]);

    std::vector<double> alpha(N);
    fftw_plan forward_plan = fftw_plan_r2r_1d(N,in.data(),alpha.data(),FFTW_REDFT10,FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);

    dataio.writeFile("fftw_alphaval.dat",alpha);

    for(auto& item : alpha)
        item /= static_cast<double>(-2.0*N);

    std::vector<double> beta(N,0.0);
    beta[N-1] = 0.0;
    beta[N-2] =  2.0*(N-1)*alpha[N-1]; 
    for(int k = N-2; k > 0; k--)
        beta[k-1] = beta[k+1] + 2.0*k*alpha[k];

    dataio.writeFile("fftw_betaval.dat",beta);
    std::vector<double> expect(N);
    fftw_plan reverse_plan = fftw_plan_r2r_1d(N,beta.data(),expect.data(),FFTW_REDFT01,FFTW_ESTIMATE);
    fftw_execute(reverse_plan);
    fftw_destroy_plan(reverse_plan);

    dataio.writeFile("fftw_der_compare.dat",in,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}





