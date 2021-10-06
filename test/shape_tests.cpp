
#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/shape/shapeT.h>
#include <pwutils/report/dataio.hpp>
#include <iostream>

TEST(SHAPE_TEST,GAUSST)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;
    std::vector<double> y(nt);

    spida::UniformGridRVT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(car_freq);
    shape.shapeRV(y);

//    pw::DataIO dataio("outfolder");
//    dataio.writeFile("gauss.dat",grid.getT(),y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));
	EXPECT_NEAR(223606797.74997896,*itmax,1e-6);
    EXPECT_NEAR(-200519096.49044815,*itmin,1e-6);
	EXPECT_NEAR(90533935.926272079,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(90533935.926272079,y[(nt+nt/16)/2],1e-6);
}

TEST(SHAPE_TEST,SECHT)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;
    std::vector<double> y(nt);

    spida::UniformGridRVT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::SechT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(car_freq);
    shape.shapeRV(y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_NEAR(223606797.74997896,*itmax,1e-6);
    EXPECT_NEAR(-211808302.73070601,*itmin,1e-6);
	EXPECT_NEAR(122726544.58500151,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(122726544.58500151,y[(nt+nt/16)/2],1e-6);
}

TEST(SHAPE_TEST,AIRYT)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -320e-15;
    double maxT = 120e-15;
    double apod = 0.2;
    std::vector<double> y(nt);

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::AiryT shape(grid,std::sqrt(I0),tp,apod);
    shape.setFastPhase(omega0);
    shape.shapeRV(y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_NEAR(111197223.83176377,*itmax,1e-6);
    EXPECT_NEAR(-115045662.1447591,*itmin,1e-6);
	EXPECT_NEAR(8518347.9235188272,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(-18030265.656660724,y[(nt+nt/16)/2],1e-6);
}

TEST(SHAPE_TEST,BESSELT)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;
    double apod = 0.5;
    std::vector<double> y(nt);

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::BesselT shape(grid,std::sqrt(I0),tp,apod);
    shape.setFastPhase(omega0);
    shape.shapeRV(y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_DOUBLE_EQ(223606797.74997896,*itmax);
    EXPECT_DOUBLE_EQ(-185314259.84787187,*itmin);
	EXPECT_DOUBLE_EQ(46643813.615054831,y[(nt-nt/16)/2]);
	EXPECT_DOUBLE_EQ(46643813.615054831,y[(nt+nt/16)/2]);
}




    







