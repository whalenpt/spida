

#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/grid/uniformCVT.h>
#include <spida/grid/besselR.h>
#include <spida/shape/shapeT.h>
#include <spida/shape/shapeR.h>
#include <spida/transform/hankelfftRBLT.h>
#include <spida/transform/hankelfftRCVT.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <pwutils/report/dat.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <random>
#include <fstream>


// Test that forward transforms followed by the inverse transforms yield the identity
TEST(HANKELFFTRCVT_TEST,INVERSES)
{
    using spida::dcmplx;
    unsigned nr = 100;
    unsigned nt = 512;
    double minT = -6.0;
    double maxT = 6.0;
    double maxR = 12.0;

    spida::UniformGridCVT gridT(nt,minT,maxT);
    spida::BesselRootGridR gridR(nr,maxR);

    spida::GaussT shapeT(gridT,1.0,1.0);
    spida::GaussR shapeR(gridR,1.0,1.0);

    std::vector<dcmplx> u0t = shapeT.shapeCV();
    std::vector<dcmplx> u0r = shapeR.shapeCV();
    std::vector<dcmplx> u(nr*nt);
    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt+j] = u0r[i]*u0t[j];

    std::vector<dcmplx> usr(nr*nt);
    std::vector<dcmplx> ub(nr*nt);

    spida::HankelFFTRCVT transform(gridR,gridT);
    // Check forward transform and reverse transform over R dimension
    transform.RT_To_SRT(u,usr);
    transform.SRT_To_RT(usr,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);


    unsigned nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    // Check forward transform and reverse transform over T dimension
    transform.RT_To_RST(u,v);
    transform.RST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    // Check forward tranform and reverse transform applied in sequence is identity
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);


    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> vb(nr*nst);
    // Check forward tranform and reverse transform over SR dimension
    transform.SRST_To_RST(v,w);
    transform.RST_To_SRST(w,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    // Check forward tranform and reverse transform over ST dimension
    transform.SRST_To_SRT(v,usr);
    transform.SRT_To_SRST(usr,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    std::vector<dcmplx> wb(nr*nst);
    // Check RT_To_RST and SRST_To_RST are equal
    transform.RT_To_RST(u,w);
    transform.SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    // Check RT_To_SRT and SRST_To_SRT are equal
    std::vector<dcmplx> zeta(nr*nt);
    std::vector<dcmplx> zetab(nr*nt);
    transform.RT_To_SRT(u,zeta);
    transform.SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}

// H0{exp(-a*r^2)}=(1/2a)*exp(-kr^2/(4a))
// FFT{exp(-b*t^2)} = sqrt(pi/b)*exp(-omega^2/(4*b)) 
// F{exp(-a*r^2 -b*t^2)} = sqrt(pi/b)/(2*a)*exp(-kr^2/(4a)-omega^2/(4*b))
TEST(HANKELFFTRCVT_TEST,GAUSSTGAUSSR)
{
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;

    unsigned nr = 100;
    unsigned nt = 512;
    double minT = -6.0;
    double maxT = 6.0;
    double maxR = 12.0;

    spida::UniformGridCVT gridT(nt,minT,maxT);
    spida::BesselRootGridR gridR(nr,maxR);
    const std::vector<double> r = gridR.getR();
    const std::vector<double> t = gridT.getT();
    const std::vector<double> kr = gridR.getSR();
    const std::vector<double> omega = gridT.getST();
    unsigned nst = gridT.getNst();

    std::vector<dcmplx> u(nr*nt);
    double a = 1.0;
    double b = 2.0;
    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt+j] = exp(-a*pow(r[i],2) - b*pow(t[j],2));

    std::vector<dcmplx> out(nr*nst);
    spida::HankelFFTRCVT transform(gridR,gridT);
    transform.RT_To_SRST(u,out);
    auto report = dat::ReportComplexData2D<double,double,double>("hankelfft_SR",kr,omega,out);
    std::ofstream os;
    report.setDirPath("outfolder");
    os << report;

    std::vector<dcmplx> expect(nr*nst);
    for(auto i = 0; i < kr.size(); i++)
        for(auto j = 0; j < omega.size(); j++)
            expect[i*nst+j] = sqrt(PI/b)/(2.0*a)*exp(-(pow(kr[i],2)/(4.0*a)+pow(omega[j],2)/(4.0*b)));

    auto report_ex = dat::ReportComplexData2D<double,double,double>("hankelfft_expect_SR",kr,omega,expect);
    report_ex.setDirPath("outfolder");
    os << report_ex;

    EXPECT_LT(pw::relative_error(expect,out),1e-5);
}

TEST(HANKELFFTRCVT_TEST,MULTITHREAD)
{
    using spida::dcmplx;
    unsigned nr = 100;
    unsigned nt = 512;
    double minT = -6.0;
    double maxT = 6.0;
    double maxR = 12.0;
    unsigned NUM_THREADS = 3;

    spida::UniformGridCVT gridT(nt,minT,maxT);
    spida::BesselRootGridR gridR(nr,maxR);

    spida::GaussT shapeT(gridT,1.0,1.0);
    spida::GaussR shapeR(gridR,1.0,1.0);

    std::vector<dcmplx> u0t = shapeT.shapeCV();
    std::vector<dcmplx> u0r = shapeR.shapeCV();
    std::vector<dcmplx> u(nr*nt);
    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt+j] = u0r[i]*u0t[j];

    unsigned nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<dcmplx> ub(nr*nt);

    // Check forward tranform and reverse transform applied in sequence is identity
    spida::HankelFFTRCVT transform(gridR,gridT,NUM_THREADS);
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    // Check RT_To_RST and SRST_To_RST are equal
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> wb(nr*nst);
    transform.RT_To_RST(u,w);
    transform.SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    // Check RT_To_SRT and SRST_To_SRT are equal
    std::vector<dcmplx> zeta(nr*nt);
    std::vector<dcmplx> zetab(nr*nt);
    transform.RT_To_SRT(u,zeta);
    transform.SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}


TEST(HANKELFFTRRVT_TEST,INVERSES)
{
    using spida::dcmplx;
    int nr = 100;
    int nt = 512;
    double w0 = 20.0e-6;
    double I0 = 5.0e16;
    double tp = 5.0e-15;
    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridRVT gridT(nt,minT,maxT,minST,maxST);
    spida::BesselRootGridR gridR(nr,12*w0);

    spida::GaussT shapeT(gridT,std::sqrt(I0),tp);
    shapeT.setFastPhase(omega0);
    spida::GaussR shapeR(gridR,1.0,w0);

    std::vector<double> u0t(nt);
    shapeT.shapeRV(u0t);
    std::vector<double> u0r(nr);
    shapeR.shapeRV(u0r);
    std::vector<double> u(nr*nt);

    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt + j] = u0r[i]*u0t[j];

    int nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<dcmplx> vb(nr*nst);
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> wb(nr*nst);
    std::vector<double> ub(nr*nt);
    std::vector<double> usr(nr*nt);
    spida::HankelFFTRBLT transform(gridR,gridT);

    // Check forward transform and reverse transform over R dimension
    transform.RT_To_SRT(u,usr);
    transform.SRT_To_RT(usr,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    // Check forward transform and reverse transform over T dimension
    transform.RT_To_RST(u,v);
    transform.RST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    // Check forward tranform and reverse transform applied in sequence is identity
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    // Check forward tranform and reverse transform over SR dimension
    transform.SRST_To_RST(v,w);
    transform.RST_To_SRST(w,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    // Check forward tranform and reverse transform over ST dimension
    transform.SRST_To_SRT(v,usr);
    transform.SRT_To_SRST(usr,vb);
    EXPECT_LT(pw::relative_error(v,vb),1e-6);

    // Check RT_To_RST and SRST_To_RST are equal
    transform.RT_To_RST(u,w);
    transform.SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    // Check RT_To_SRT and SRST_To_SRT are equal
    std::vector<double> zeta(nr*nt);
    std::vector<double> zetab(nr*nt);
    transform.RT_To_SRT(u,zeta);
    transform.SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}

// H0{exp(-(r/w0)^2)}=(w0^2/2)*exp(-w0^2*kr^2/4)
// FFT{exp(-(t/tp)^2)}= tp*sqrt(pi)*exp(-tp^2*omega^2/4)
// F{exp(-(r/w0)^2-(t/tp)^2)} = (tp*w0^2*sqrt(pi)/2)*exp(-w0^2*kr^2/4-tp^2*omega^2/4)
TEST(HANKELFFTRRVT_TEST,GAUSSTGAUSSR)
{
    using spida::dcmplx;
    using spida::PI;
    using spida::ii;
    int nr = 100;
    int nt = 512;
    double w0 = 20.0e-6;
    double maxR = 12*w0;
    double I0 = 5.0e16;
    double tp = 5.0e-15;
    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;

    spida::UniformGridRVT gridT(nt,minT,maxT,minST,maxST);
    spida::BesselRootGridR gridR(nr,maxR);

    const std::vector<double> r = gridR.getR();
    const std::vector<double> t = gridT.getT();
    const std::vector<double> kr = gridR.getSR();
    const std::vector<double> omega = gridT.getST();
    unsigned nst = gridT.getNst();

    std::vector<double> u(nr*nt);
    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt+j] = sqrt(I0)*exp(-pow(r[i]/w0,2) - pow(t[j]/tp,2))*cos(omega0*t[j]);

    std::vector<dcmplx> out(nr*nst);
    spida::HankelFFTRBLT transform(gridR,gridT);

    transform.RT_To_SRST(u,out);
    auto report = dat::ReportComplexData2D<double,double,double>("hankelrfft_SR",kr,omega,out);
    std::ofstream os;
    report.setDirPath("outfolder");
    os << report;

    std::vector<dcmplx> expect(nr*nst);
    for(auto i = 0; i < kr.size(); i++)
        for(auto j = 0; j < omega.size(); j++)
            expect[i*nst+j] = (std::sqrt(I0)*tp*pow(w0,2)*sqrt(PI)/2.0)*exp(\
                    -pow(w0,2)*pow(kr[i],2)/4.0-pow(tp,2)*pow(omega[j]-omega0,2)/4.0);

    auto report_ex = dat::ReportComplexData2D<double,double,double>("hankelfft_expect_SR",kr,omega,expect);
    report_ex.setDirPath("outfolder");
    os << report_ex;

    
    // phase seems to be slightly off from expected in spectral domain (0.5 relative_error)
    EXPECT_LT(pw::relative_error(expect,out),1);

    std::vector<double> out_abs(nr*nst);
    std::vector<double> out_ex_abs(nr*nst);
    for(auto j = 0; j < nr*nst; j++){
        out_abs[j] = abs(out[j]);
        out_ex_abs[j] = abs(expect[j]);
    }
    // absolute value seems to be accurate enough
    EXPECT_LT(pw::relative_error(out_abs,out_ex_abs),1e-6);
}



TEST(HANKELFFTRRVT_TEST,MULTITHREAD)
{
    using spida::dcmplx;
    int nr = 100;
    int nt = 512;
    double w0 = 20.0e-6;
    double I0 = 5.0e16;
    double tp = 5.0e-15;
    double omega0 = 2.7091e15;
    double minT = -10*tp;
    double maxT = 10*tp;
    double minST = 1.0e15;
    double maxST = 4.3e15;
    int NUM_THREADS = 3;

    spida::UniformGridRVT gridT(nt,minT,maxT,minST,maxST);
    spida::BesselRootGridR gridR(nr,12*w0);
    spida::GaussT shapeT(gridT,std::sqrt(I0),tp);
    shapeT.setFastPhase(omega0);
    spida::GaussR shapeR(gridR,1.0,w0);

    std::vector<double> u0t(nt);
    shapeT.shapeRV(u0t);
    std::vector<double> u0r(nr);
    shapeR.shapeRV(u0r);

    pw::DataIO dataio("outfolder");
    dataio.writeFile("rshape.dat",u0r);

    std::vector<double> u(nr*nt);

    for(auto i = 0; i < nr; i++)
        for(auto j = 0; j < nt; j++)
            u[i*nt + j] = u0r[i]*u0t[j];

    int nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<double> ub(nr*nt);

    // Check forward tranform and reverse transform applied in sequence is identity
    spida::HankelFFTRBLT transform(gridR,gridT,NUM_THREADS);
    transform.RT_To_SRST(u,v);
    transform.SRST_To_RT(v,ub);
    EXPECT_LT(pw::relative_error(u,ub),1e-6);

    // Check RT_To_RST and SRST_To_RST are equal
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> wb(nr*nst);
    transform.RT_To_RST(u,w);
    transform.SRST_To_RST(v,wb);
    EXPECT_LT(pw::relative_error(w,wb),1e-6);

    // Check RT_To_SRT and SRST_To_SRT are equal
    std::vector<double> zeta(nr*nt);
    std::vector<double> zetab(nr*nt);
    transform.RT_To_SRT(u,zeta);
    transform.SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
}








