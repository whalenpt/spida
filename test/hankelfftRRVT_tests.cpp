#include <algorithm>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>
#include <gtest/gtest.h>
#include <pwutils/pwmath.hpp>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/besselR.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelfftRRVT.h>


class HankelFFT_RRVTTEST : public ::testing::Test 
{
    protected:
        void SetUp() override {
            nt = 512;
            tp = 5.0e-15;
            minT = -10*tp; maxT = 10*tp;
            minST = 1.0e15; maxST = 4.3e15;
            nr = 100;
            w0 = 20.0e-6;
            maxR = 12.0*w0;
            I0 = 5.0e16;
            omega0 = 2.7091e15;

            threads = 2;

            // setup grids and transform
            gridT = std::make_unique<spida::UniformGridRVT>(nt,minT,maxT,minST,maxST);
            nst = gridT->getNst();
            gridR = std::make_unique<spida::BesselRootGridR>(nr,maxR);
            transform = std::make_unique<spida::HankelFFTRRVT>(*gridR,*gridT);
            transform_threaded = std::make_unique<spida::HankelFFTRRVT>(*gridR,*gridT,threads);

            u.resize(nr*nt);
            v.resize(nr*nst);
            initGaussRT(u);
        }

        // initialize arrays as gauss-gauss field
        void initGaussRT(std::vector<double>& arr) const{

            spida::GaussT shapeT(*gridT,std::sqrt(I0),tp);
            shapeT.setFastPhase(omega0);
            spida::GaussR shapeR(*gridR,1.0,w0);

            auto u0t = shapeT.shapeRV();
            auto u0r = shapeR.shapeRV();
            for(size_t i = 0; i < nr; i++)
                for(size_t j = 0; j < nt; j++)
                    arr[i*nt + j] = u0r[i]*u0t[j];
        }

        unsigned nt; 
        unsigned nr;
        unsigned nst;
        double minT; 
        double maxT; 
        double minST; 
        double maxST;
        double maxR;
        double tp; 
        double w0; 
        double I0; 
        double omega0;
        unsigned threads;
        std::unique_ptr<spida::HankelFFTRRVT> transform;
        std::unique_ptr<spida::HankelFFTRRVT> transform_threaded;
        std::unique_ptr<spida::UniformGridRVT> gridT;
        std::unique_ptr<spida::BesselRootGridR> gridR;
        std::vector<double> u;
        std::vector<spida::dcmplx> v;
};


TEST_F(HankelFFT_RRVTTEST,INVERSES1)
{
    using spida::dcmplx;
    std::vector<double> ub(nr*nt);
    transform->RT_To_SRST(u,v);
    transform->SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,INVERSES2)
{
    using spida::dcmplx;
    std::vector<double> usr(nr*nt);
    std::vector<double> ub(nr*nt);

    transform->RT_To_SRT(u,usr);
    transform->SRT_To_RT(usr,ub); 
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    transform->RT_To_RST(u,v);
    transform->RST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,INVERSES3)
{
    using spida::dcmplx;
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> vb(nr*nst);
    std::vector<double> usr(nr*nt);

    transform->RT_To_SRST(u,v);
    transform->SRST_To_RST(v,w);
    transform->RST_To_SRST(w,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    transform->SRST_To_SRT(v,usr);
    transform->SRT_To_SRST(usr,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,INVERSES4)
{
    using spida::dcmplx;
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> wb(nr*nst);
    std::vector<double> zeta(nr*nt);
    std::vector<double> zetab(nr*nt);

    transform->RT_To_SRST(u,v);
    transform->RT_To_RST(u,w);
    transform->SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    transform->RT_To_SRT(u,zeta);
    transform->SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,MULTITHREADED1)
{
    using spida::dcmplx;
    std::vector<double> ub(nr*nt);
    transform_threaded->RT_To_SRST(u,v);
    transform_threaded->SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,MULTITHREADED2)
{
    using spida::dcmplx;
    std::vector<double> usr(nr*nt);
    std::vector<double> ub(nr*nt);

    transform_threaded->RT_To_SRT(u,usr);
    transform_threaded->SRT_To_RT(usr,ub); 
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    transform_threaded->RT_To_RST(u,v);
    transform_threaded->RST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,MULTITHREADED3)
{
    using spida::dcmplx;
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> vb(nr*nst);
    std::vector<double> usr(nr*nt);

    transform_threaded->RT_To_SRST(u,v);
    transform_threaded->SRST_To_RST(v,w);
    transform_threaded->RST_To_SRST(w,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    transform_threaded->SRST_To_SRT(v,usr);
    transform_threaded->SRT_To_SRST(usr,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);
}

TEST_F(HankelFFT_RRVTTEST,MULTITHREADED4)
{
    using spida::dcmplx;
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> wb(nr*nst);
    std::vector<double> zeta(nr*nt);
    std::vector<double> zetab(nr*nt);

    transform_threaded->RT_To_SRST(u,v);
    transform_threaded->RT_To_RST(u,w);
    transform_threaded->SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    transform_threaded->RT_To_SRT(u,zeta);
    transform_threaded->SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}

// H0{exp(-(r/w0)^2)}=(w0^2/2)*exp(-w0^2*kr^2/4)
// FFT{exp(-(t/tp)^2)}= tp*sqrt(pi)*exp(-tp^2*omega^2/4)
// F{exp(-(r/w0)^2-(t/tp)^2)} = (tp*w0^2*sqrt(pi)/2)*exp(-w0^2*kr^2/4-tp^2*omega^2/4)
TEST_F(HankelFFT_RRVTTEST,GAUSSTGAUSSR)
{
    using spida::dcmplx;
    using spida::PI;
    std::vector<double> in(nr*nt);
    std::vector<dcmplx> out(nr*nst);
    const auto& r = gridR->getR();
    const auto& t = gridT->getT();
    const auto& kr = gridR->getSR();
    const auto& omega = gridT->getST();

    for(unsigned i = 0; i < nr; i++)
        for(unsigned j = 0; j < nt; j++)
            in[i*nt+j] = sqrt(I0)*exp(-pow(r[i]/w0,2) - pow(t[j]/tp,2))*cos(omega0*t[j]);

    transform->RT_To_SRST(in,out);

    // y = f(t)*cos(i\omega0t) - > FFT{y} = (FFT{f(\omega - \omega0)}+FFT{f(\omega+\omega0)})/2
    // For real fields, fft taken over positive frequencies: FFT_real{y} = FFT_real{f(\omega-\omega0)}/2  
    std::vector<dcmplx> expect(nr*nst);
    for(size_t i = 0; i < kr.size(); i++)
        for(size_t j = 0; j < omega.size(); j++)
            expect[i*nst+j] = 0.5*(std::sqrt(I0)*tp*pow(w0,2)*sqrt(PI)/2.0)*exp(\
                    -pow(w0,2)*pow(kr[i],2)/4.0-pow(tp,2)*pow(omega[j]-omega0,2)/4.0);

    EXPECT_LT(pw::relative_error(expect,out),1.0e-5);
}
