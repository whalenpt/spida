

#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformCVT.h>
#include <spida/grid/besselR.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelfftRRVT.h>
#include <spida/transform/hankelfftRCVT.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <pwutils/report/dat.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <random>
#include <fstream>

class HankelFFT_RCVTTEST : public ::testing::Test 
{
    protected:
        void SetUp() override {
            nt = 512;
            minT = -6.0; maxT = 6.0;
            nr = 100;
            maxR = 12.0;
            threads = 2;

            // setup grids and transform
            gridT = new spida::UniformGridCVT(nt,minT,maxT);
            nst = gridT->getNst();
            gridR = new spida::BesselRootGridR(nr,maxR);
            transform = new spida::HankelFFTRCVT(*gridR,*gridT);
            transform_threaded = new spida::HankelFFTRCVT(*gridR,*gridT,threads);
            u.resize(nr*nt);
            v.resize(nr*nst);
            initGaussRT(u);
        }

        void TearDown() override {
            delete transform;
            delete transform_threaded;
            delete gridT;
            delete gridR;
        }

        // initial field
        void initGaussRT(std::vector<spida::dcmplx>& arr){
            spida::GaussT shapeT(*gridT,1.0,1.0);
            spida::GaussR shapeR(*gridR,1.0,1.0);
            std::vector<spida::dcmplx> u0t = shapeT.shapeCV();
            std::vector<spida::dcmplx> u0r = shapeR.shapeCV();
            for(auto i = 0; i < nr; i++)
                for(auto j = 0; j < nt; j++)
                    arr[i*nt+j] = u0r[i]*u0t[j];
        }

        unsigned nr, nt, nst;
        double minT, maxT;
        double maxR;
        unsigned threads;
        spida::HankelFFTRCVT* transform;
        spida::HankelFFTRCVT* transform_threaded;
        spida::UniformGridCVT* gridT;
        spida::BesselRootGridR* gridR;
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
    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            in[i*nt+j] = exp(-a*pow(r[i],2) - b*pow(t[j],2));

    transform->RT_To_SRST(in,out);
    auto report = dat::ReportComplexData2D<double,double,double>(\
        "hankelfft_SR",kr,omega,out);
    std::ofstream os;
    report.setDirPath("outfolder");
    os << report;

    std::vector<dcmplx> expect(nr*nst);
    for(auto i = 0; i < kr.size(); i++){
        for(auto j = 0; j < omega.size(); j++){
            expect[i*nst+j] = sqrt(spida::PI/b)/(2.0*a)*exp(\
            -(pow(kr[i],2)/(4.0*a)+pow(omega[j],2)/(4.0*b)));
        }
    }

    auto report_ex = dat::ReportComplexData2D<double,double,double>(\
        "hankelfft_expect_SR",kr,omega,expect);
    report_ex.setDirPath("outfolder");
    os << report_ex;

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
            gridT = new spida::UniformGridRVT(nt,minT,maxT,minST,maxST);
            nst = gridT->getNst();
            gridR = new spida::BesselRootGridR(nr,maxR);
            transform = new spida::HankelFFTRRVT(*gridR,*gridT);
            transform_threaded = new spida::HankelFFTRRVT(*gridR,*gridT,threads);

            u.resize(nr*nt);
            v.resize(nr*nst);
            initGaussRT(u);
        }


        void TearDown() override {
            delete transform;
            delete transform_threaded;
            delete gridT;
            delete gridR;
        }

        // initialize arrays as gauss-gauss field
        void initGaussRT(std::vector<double>& arr){

            spida::GaussT shapeT(*gridT,std::sqrt(I0),tp);
            shapeT.setFastPhase(omega0);
            spida::GaussR shapeR(*gridR,1.0,w0);

            std::vector<double> u0t(nt);
            shapeT.shapeRV(u0t);
            std::vector<double> u0r(nr);
            shapeR.shapeRV(u0r);
            for(size_t i = 0; i < nr; i++)
                for(size_t j = 0; j < nt; j++)
                    arr[i*nt + j] = u0r[i]*u0t[j];
        }

        unsigned nt, nr, nst;
        double minT, maxT, minST, maxST;
        double maxR;
        double tp, w0, I0, omega0;
        unsigned threads;
        spida::HankelFFTRRVT* transform;
        spida::HankelFFTRRVT* transform_threaded;
        spida::UniformGridRVT* gridT;
        spida::BesselRootGridR* gridR;
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
    const std::vector<double>& r = gridR->getR();
    const std::vector<double>& t = gridT->getT();
    const std::vector<double>& kr = gridR->getSR();
    const std::vector<double>& omega = gridT->getST();

    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            in[i*nt+j] = sqrt(I0)*exp(-pow(r[i]/w0,2) - pow(t[j]/tp,2))*cos(omega0*t[j]);

    transform->RT_To_SRST(in,out);
    auto report = dat::ReportComplexData2D<double,double,double>("hankelrfft_SR",kr,omega,out);
    std::ofstream os;
    report.setDirPath("outfolder");
    os << report;

    // y = f(t)*cos(i\omega0t) - > FFT{y} = (FFT{f(\omega - \omega0)}+FFT{f(\omega+\omega0)})/2
    // For real fields, fft taken over positive frequencies: FFT_real{y} = FFT_real{f(\omega-\omega0)}/2  
    std::vector<dcmplx> expect(nr*nst);
    for(auto i = 0; i < kr.size(); i++)
        for(auto j = 0; j < omega.size(); j++)
            expect[i*nst+j] = 0.5*(std::sqrt(I0)*tp*pow(w0,2)*sqrt(PI)/2.0)*exp(\
                    -pow(w0,2)*pow(kr[i],2)/4.0-pow(tp,2)*pow(omega[j]-omega0,2)/4.0);

    auto report_ex = dat::ReportComplexData2D<double,double,double>(\
        "hankelrfft_expect_SR",kr,omega,expect);
    report_ex.setDirPath("outfolder");
    os << report_ex;

    EXPECT_LT(pw::relative_error(expect,out),1.0e-5);
}
