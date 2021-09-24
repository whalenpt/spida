
#include <gtest/gtest.h>
#include <spida/grid/uniformX.h>
#include <spida/shape/shapeX.h>
#include <spida/transform/fftX.h>
#include <spida/SpidaX.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <random>

TEST(FFT_TEST,INVERSES)
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

    spida::FFTX tr(spida::UniformGridX{N,-1,1});
    tr.X_To_SX(in,out);
    tr.SX_To_X(out,expect);

    dataio.writeFile("fft_check.dat",in,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

TEST(FFT_TEST,GAUSS)
{
	unsigned N = 256;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    spida::UniformGridX grid(N,-6,6);
    const std::vector<double> x = grid.getX();
    const std::vector<double> kx = grid.getSX();
    double alpha = 1.0;
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(-alpha*pow(x[i],2));
    for(auto i = 0; i < kx.size(); i++)
        expect[i] = sqrt(PI/alpha)*exp(-pow(PI*kx[i],2)/alpha);
    dataio.writeFile("fft_gauss.dat",in);

    spida::FFTX tr(grid);
    tr.X_To_SX(in,out);

    std::vector<double> reals_out(N);
    for(auto i = 0; i < out.size(); i++)
        reals_out[i] = out[i].real();
    std::vector<double> reals_expect(N);
    for(auto i = 0; i < expect.size(); i++)
        reals_expect[i] = expect[i].real();

    dataio.writeFile("fft_gauss_check.dat",reals_expect,reals_out);
    EXPECT_LT(pw::relative_error(out,expect),1e-6);
}

TEST(FFT_TEST,COS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    spida::UniformGridX grid(N,0,2*spida::PI);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = cos(8*x[i]);

    spida::FFTX tr(grid);
    tr.X_To_SX(in,out);
    dataio.writeFile("fft_cos_check.dat",out);

	EXPECT_DOUBLE_EQ(out[8].real(),static_cast<double>(N)/2);
}



TEST(FFT_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");

    using spida::dcmplx;
    using spida::PI;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridX grid(N,0,2*PI);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = sin(x[i]);
    for(auto i = 0; i < x.size(); i++)
        expect[i] = cos(x[i]);

    spida::SpidaX spidaX{grid};
    spidaX.dX(in,out);
    dataio.writeFile("fft_der_sin.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

TEST(FFT_TEST,DERIVATIVE_GAUSS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");

    using spida::dcmplx;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridX grid(N,-6,6);
    const std::vector<double> x = grid.getX();
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(-pow(x[i],2));
    for(auto i = 0; i < x.size(); i++)
        expect[i] = -2.0*x[i]*exp(-pow(x[i],2));


    spida::SpidaX spidaX{grid};
    spidaX.dX(in,out);
    dataio.writeFile("fft_der_gauss.dat",expect,out);
    dataio.writeFile("fft_x_sx.dat",grid.getX(),grid.getSX());
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}



/*
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
*/





