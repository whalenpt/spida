#include <algorithm>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>
#include <gtest/gtest.h>
#include <pwutils/pwmath.hpp>
#include <spida/helper/constants.h>
#include <spida/grid/uniformCVT.h>
#include <spida/grid/besselR.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelfftRCVT.h>

class HankelFFT_RCVTTEST : public ::testing::Test 
{
    protected:
        void SetUp() override {
            nt = 512;
            minT = -6.0; 
            maxT = 6.0;
            nr = 100;
            maxR = 12.0;
            threads = 2;

            // setup grids and transform
            gridT = std::make_unique<spida::UniformGridCVT>(nt,minT,maxT);
            nst = gridT->getNst();
            gridR = std::make_unique<spida::BesselRootGridR>(nr,maxR);
            transform = std::make_unique<spida::HankelFFTRCVT>(*gridR,*gridT);
            transform_threaded = std::make_unique<spida::HankelFFTRCVT>(*gridR,*gridT,threads);
            u.resize(nr*nt);
            v.resize(nr*nst);
            initGaussRT(u);
        }

        void TearDown() override {}

        // initial field
        void initGaussRT(std::vector<spida::dcmplx>& arr) const{
            spida::GaussT shapeT(*gridT,1.0,1.0);
            spida::GaussR shapeR(*gridR,1.0,1.0);
            std::vector<spida::dcmplx> u0t = shapeT.shapeCV();
            std::vector<spida::dcmplx> u0r = shapeR.shapeCV();
            for(unsigned i = 0; i < nr; i++)
                for(unsigned j = 0; j < nt; j++)
                    arr[i*nt+j] = u0r[i]*u0t[j];
        }

        unsigned nr;
        unsigned nt;
        unsigned nst;
        double minT; 
        double maxT;
        double maxR;
        unsigned threads;
        std::unique_ptr<spida::HankelFFTRCVT> transform;
        std::unique_ptr<spida::HankelFFTRCVT> transform_threaded;
        std::unique_ptr<spida::UniformGridCVT> gridT;
        std::unique_ptr<spida::BesselRootGridR> gridR;
        std::vector<spida::dcmplx> u;
        std::vector<spida::dcmplx> v;
};

// Test RT_To_SRST and SRST_To_RT transforms are inverses
TEST_F(HankelFFT_RCVTTEST,INVERSES)
{
    using spida::dcmplx;
    std::vector<dcmplx> ub(nr*nt);
    transform->RT_To_SRST(u,v);
    transform->SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

// Check RT_To_SRT and SRT_To_RT transforms are inverses
// Check RT_To_RST and RST_To_RT transforms are inverses
TEST_F(HankelFFT_RCVTTEST,INVERSES2)
{
    std::vector<spida::dcmplx> usr(nr*nt);
    std::vector<spida::dcmplx> ub(nr*nt);
    // Check forward transform and reverse transform over R dimension
    transform->RT_To_SRT(u,usr);
    transform->SRT_To_RT(usr,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
    // Check forward transform and reverse transform over T dimension
    transform->RT_To_RST(u,v);
    transform->RST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

// Check SRST_To_RST and RST_To_SRST transforms are inverses
// Check SRST_To_SRT and SRT_To_SRST transforms are inverses
TEST_F(HankelFFT_RCVTTEST,INVERSES3)
{
    std::vector<spida::dcmplx> w(nr*nst);
    std::vector<spida::dcmplx> vb(nr*nst);
    std::vector<spida::dcmplx> usr(nr*nt);

    transform->RT_To_SRST(u,v);

    // Check forward tranform and reverse transform over SR dimension
    transform->SRST_To_RST(v,w);
    transform->RST_To_SRST(w,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    // Check forward tranform and reverse transform over ST dimension
    transform->SRST_To_SRT(v,usr);
    transform->SRT_To_SRST(usr,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);
}

// Check RT_To_RST and SRST_To_RST are equal
// Check RT_To_SRT and SRST_To_SRT are equal
TEST_F(HankelFFT_RCVTTEST,INVERSES4)
{
    std::vector<spida::dcmplx> w(nr*nst);
    std::vector<spida::dcmplx> wb(nr*nst);
    transform->RT_To_SRST(u,v);

    // Check RT_To_RST and SRST_To_RST are equal
    transform->RT_To_RST(u,w);
    transform->SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    // Check RT_To_SRT and SRST_To_SRT are equal
    std::vector<spida::dcmplx> zeta(nr*nt);
    std::vector<spida::dcmplx> zetab(nr*nt);
    transform->RT_To_SRT(u,zeta);
    transform->SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}


// H0{exp(-a*r^2)}=(1/2a)*exp(-kr^2/(4a))
// FFT{exp(-b*t^2)} = sqrt(pi/b)*exp(-omega^2/(4*b)) 
// F{exp(-a*r^2 -b*t^2)} = sqrt(pi/b)/(2*a)*exp(-kr^2/(4a)-omega^2/(4*b))
TEST_F(HankelFFT_RCVTTEST,GAUSSTGAUSSR)
{
    using spida::dcmplx;
    std::vector<dcmplx> in(nr*nt);
    std::vector<dcmplx> out(nr*nst);
    double a = 1.0;
    double b = 2.0;
    const std::vector<double>& r = gridR->getR();
    const std::vector<double>& t = gridT->getT();
    const std::vector<double>& kr = gridR->getSR();
    const std::vector<double>& omega = gridT->getST();
    for(unsigned i = 0; i < nr; i++)
        for(unsigned j = 0; j < nt; j++)
            in[i*nt+j] = exp(-a*pow(r[i],2) - b*pow(t[j],2));

    transform->RT_To_SRST(in,out);

    std::vector<dcmplx> expect(nr*nst);
    for(size_t i = 0; i < kr.size(); i++){
        for(size_t j = 0; j < omega.size(); j++){
            expect[i*nst+j] = sqrt(spida::PI/b)/(2.0*a)*exp(\
            -(pow(kr[i],2)/(4.0*a)+pow(omega[j],2)/(4.0*b)));
        }
    }

    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

// Test RT_To_SRST and SRST_To_RT transforms are inverses (multithread)
TEST_F(HankelFFT_RCVTTEST,MULTITHREAD1)
{
    using spida::dcmplx;
    std::vector<dcmplx> ub(nr*nt);
    // Check forward tranform and reverse transform applied in sequence is identity
    transform_threaded->RT_To_SRST(u,v);
    transform_threaded->SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

// Check RT_To_SRT and SRT_To_RT transforms are inverses (multithread)
// Check RT_To_RST and RST_To_RT transforms are inverses (multithread)
TEST_F(HankelFFT_RCVTTEST,MULTITHREAD2)
{
    std::vector<spida::dcmplx> usr(nr*nt);
    std::vector<spida::dcmplx> ub(nr*nt);
    // Check forward transform and reverse transform over R dimension
    transform_threaded->RT_To_SRT(u,usr);
    transform_threaded->SRT_To_RT(usr,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
    // Check forward transform and reverse transform over T dimension
    transform_threaded->RT_To_RST(u,v);
    transform_threaded->RST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);
}

// Check SRST_To_RST and RST_To_SRST transforms are inverses (multithread)
// Check SRST_To_SRT and SRT_To_SRST transforms are inverses (multithread)
TEST_F(HankelFFT_RCVTTEST,MULTITHREAD3)
{
    std::vector<spida::dcmplx> w(nr*nst);
    std::vector<spida::dcmplx> vb(nr*nst);
    std::vector<spida::dcmplx> usr(nr*nt);
    transform_threaded->RT_To_SRST(u,v);
    // Check forward tranform and reverse transform over SR dimension
    transform_threaded->SRST_To_RST(v,w);
    transform_threaded->RST_To_SRST(w,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    // Check forward tranform and reverse transform over ST dimension
    transform_threaded->SRST_To_SRT(v,usr);
    transform_threaded->SRT_To_SRST(usr,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);
}