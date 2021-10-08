

#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformCVT.h>
#include <spida/shape/shapeT.h>
#include <spida/transform/fftBLT.h>
#include <spida/transform/fftCVT.h>
#include <spida/SpidaCVT.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <random>


// FFTCVT defined such that F{f(t)} = \integral_{-\inf}^{\inf}f(t)exp(i*omega*t) dt
// Test that forward fft followed by inverse fft yields identity
TEST(FFTCVT_TEST,INVERSES)
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

    spida::FFTCVT tr(spida::UniformGridCVT{N,-1,1});
    tr.T_To_ST(in,out);
    tr.ST_To_T(out,expect);

    dataio.writeFile("fftcvt_check.dat",in,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

// a = 1/tp^2
// F{exp(-a*t^2)} = sqrt(pi/a)*exp(-omega^2/(4*a)) 
// F{exp(-(t/tp)^2)}= tp*sqrt(pi)*exp(-tp^2*omega^2/4)
TEST(FFTCVT_TEST,GAUSS)
{
	unsigned N = 64;
    pw::DataIO dataio("outfolder");
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
    for(auto i = 0; i < t.size(); i++)
        in[i] = exp(-a*pow(t[i],2));
    for(auto i = 0; i < omega.size(); i++)
        expect[i] = sqrt(PI/a)*exp(-pow(omega[i],2)/(4.0*a));

    spida::FFTCVT tr(grid);
    tr.T_To_ST(in,out);
    dataio.writeFile("fftcvt_gauss.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

// F{cos(at)} = PI*(delta(omega-a) + delta(omega+a))
TEST(FFTCVT_TEST,COS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    spida::UniformGridCVT grid(N,0.0,2.0*PI);
    const std::vector<double> t = grid.getT();
    for(auto i = 0; i < t.size(); i++)
        in[i] = cos(8*t[i]);

    spida::FFTCVT tr(grid);
    tr.T_To_ST(in,out);
    dataio.writeFile("fftcvt_cos_check.dat",out);
	EXPECT_DOUBLE_EQ(out[8].real(),PI);
}

TEST(FFTCVT_TEST,DERIVATIVE_SIN)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");

    using spida::dcmplx;
    using spida::PI;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVT grid(N,0,2*PI);
    const std::vector<double> t = grid.getT();
    for(auto i = 0; i < t.size(); i++)
        in[i] = sin(t[i]);
    for(auto i = 0; i < t.size(); i++)
        expect[i] = cos(t[i]);

    spida::SpidaCVT spi{grid};
    spi.dT(in,out);
    dataio.writeFile("fftcv_der_sin.dat",expect,out);
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

TEST(FFTCVT_TEST,DERIVATIVE_GAUSS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");

    using spida::dcmplx;
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);

    spida::UniformGridCVT grid(N,-6,6);
    const std::vector<double> t = grid.getT();
    for(auto i = 0; i < t.size(); i++)
        in[i] = exp(-pow(t[i],2));
    for(auto i = 0; i < t.size(); i++)
        expect[i] = -2.0*t[i]*exp(-pow(t[i],2));

    spida::SpidaCVT spi{grid};
    spi.dT(in,out);
    dataio.writeFile("fftcv_der_gauss.dat",expect,out);
    dataio.writeFile("fftcv_t_st.dat",grid.getT(),grid.getST());
    EXPECT_LT(pw::relative_error(expect,out),1e-6);
}

/*

TEST(FFTRVT_TEST,GAUSST)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    std::vector<double> y(nt);
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());

    spida::FFTBLT transform(grid);
    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);
    shape.shapeRV(y);

    transform.T_To_ST(y,ysp);
    transform.ST_To_T(ysp,yinv);

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),33820027093.103012);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}

TEST(FFTRVT_TEST,GAUSST_POINTERS)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::FFTBLT transform(grid);

    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);

    std::vector<double> y(nt);
    std::vector<double> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());

    shape.ShapeT::shapeRV(y);
    transform.T_To_ST(y.data(),ysp.data());
    transform.ST_To_T(ysp.data(),yinv.data());

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),33820027093.103012);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}


TEST(FFTRVT_TEST,COMPLEX_GAUSST)
{
    using spida::dcmplx;
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::FFTBLT transform(grid);

    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(omega0);

    std::vector<dcmplx> y(nt);
    std::vector<dcmplx> yinv(nt);
    std::vector<dcmplx> ysp(grid.getNst());
    shape.shapeCV(y);

    transform.CVT_To_ST(y,ysp);
    transform.ST_To_CVT(ysp,yinv);

    auto maxval = pw::max(ysp);
    auto maxpos = pw::argmax(ysp);
	EXPECT_DOUBLE_EQ(abs(maxval),2*33820027093.103012);
	EXPECT_EQ(maxpos,28);
    EXPECT_LT(pw::relative_error(y,yinv),1e-6);
}
*/




/*
TEST(FFTCVT_TEST,GAUSS)
{
	unsigned N = 32;
    pw::DataIO dataio("outfolder");
    using spida::dcmplx;
    using spida::PI;

    std::vector<dcmplx> in(N,0.0);
    std::vector<dcmplx> out(N,0.0);
    std::vector<dcmplx> expect(N,0.0);

    spida::UniformGridX grid(N,-6,6);
    double L = grid.getLX();
    const std::vector<double> x = grid.getX();
    const std::vector<double> kx = grid.getSX();
    
    double alpha = 2.0;
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(-alpha*pow(x[i],2));
    for(auto i = 0; i < kx.size(); i++)
        expect[i] = (N/L)*sqrt(PI/alpha)*exp(-pow(kx[i],2)/(4*alpha));

        //expect[i] = (N/L)*sqrt(PI/alpha)*exp(-2.88202*pow(PI*L*kx[i]/static_cast<double>(N),2)/alpha);


    spida::FFTCVX tr(grid);
    tr.X_To_SX(in,out);

    dataio.writeFile("fft_gauss.dat",expect,out);

    EXPECT_LT(pw::relative_error(out,expect),1e-5);

    std::vector<double> reals_out(N);
    for(auto i = 0; i < out.size(); i++)
        reals_out[i] = out[i].real();
    std::vector<double> reals_expect(N);
    for(auto i = 0; i < expect.size(); i++)
        reals_expect[i] = expect[i].real();

    std::vector<dcmplx> outfftw(N,0.0);
    fftw_plan forward_plan = fftw_plan_dft_1d(in.size(),\
            reinterpret_cast<fftw_complex*>(in.data()),\
            reinterpret_cast<fftw_complex*>(outfftw.data()),\
            FFTW_FORWARD,\
            FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);

//    backward_plan = fftw_plan_dft_1d(in.size(),reinterpret_cast<fftw_complex*>(out.data()),\

    std::vector<double> fftw_real(N);
    for(auto i = 0; i < outfftw.size(); i++)
        fftw_real[i] = outfftw[i].real();

    dataio.writeFile("fft_gauss_check.dat",reals_expect,reals_out);
    dataio.writeFile("fft_gauss_check_fftw.dat",reals_expect,fftw_real);
    dataio.writeFile("fft_gauss_check_dcmplx.dat",expect,out);
    dataio.writeFile("fft_gauss_check_fftw_dcmplx.dat",expect,outfftw);

    std::vector<dcmplx> in2(N,0.0);
    std::vector<double> orderedx(N,0.0);
    //for(auto i = 0; i < x.size(); i++)
    for(auto i = 0; i < x.size()/2; i++)
        orderedx[i] = x[i+x.size()/2];
    for(auto i = x.size()/2; i< x.size(); i++)
        orderedx[i] = x[i-x.size()/2];
    for(auto i = 0; i < orderedx.size(); i++)
        in2[i] = exp(-alpha*pow(orderedx[i],2));
    dataio.writeFile("fft_gauss_ordered_x.dat",orderedx,x);

    std::vector<dcmplx> outfftw2(N,0.0);
    fftw_plan forward_plan2 = fftw_plan_dft_1d(in2.size(),\
            reinterpret_cast<fftw_complex*>(in2.data()),\
            reinterpret_cast<fftw_complex*>(outfftw2.data()),\
            FFTW_FORWARD,\
            FFTW_ESTIMATE);
    fftw_execute(forward_plan2);
    fftw_destroy_plan(forward_plan2);
    dataio.writeFile("fft_gauss_check_ordered.dat",expect,outfftw2);
    EXPECT_LT(pw::relative_error(outfftw2,expect),1e-5);
}
*/



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









