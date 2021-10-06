

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
#include <algorithm>
#include <numeric>
#include <functional>
#include <random>

TEST(HANKELFFTRCVT_TEST,GAUSSTGAUSSR)
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

    unsigned nst = gridT.getNst();
    std::vector<dcmplx> v(nr*nst);
    std::vector<dcmplx> vb(nr*nst);
    std::vector<dcmplx> w(nr*nst);
    std::vector<dcmplx> wb(nr*nst);
    std::vector<dcmplx> ub(nr*nt);
    std::vector<dcmplx> usr(nr*nt);
    // Check forward transform and reverse transform over R dimension

    spida::HankelFFTRCVT transform(gridR,gridT);
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
    std::vector<dcmplx> zeta(nr*nt);
    std::vector<dcmplx> zetab(nr*nt);
    transform.RT_To_SRT(u,zeta);
    transform.SRST_To_SRT(v,zetab);
    EXPECT_LT(pw::relative_error(zeta,zetab),1e-6);
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


TEST(HANKELFFTRRVT_TEST,GAUSSTGAUSSR)
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








