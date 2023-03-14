
#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformRVT.h>
#include <spida/shape/shapeT.h>
#include <iostream>


TEST(SHAPE_TEST, GAUSST_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    spida::GaussT shape(grid, A, width);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
}

TEST(SHAPE_TEST, SECHT_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    spida::SechT shape(grid, A, width);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
}

TEST(SHAPE_TEST, AIRYT_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    double apod = 1.0;
    spida::AiryT shape(grid, A, width, apod);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
    EXPECT_DOUBLE_EQ(apod, shape.apodization());
}

TEST(SHAPE_TEST, SUPERGAUSST_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    double M = 2.0;
    spida::SuperGaussT shape(grid, A, width, M);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
    EXPECT_DOUBLE_EQ(M, shape.M());
}


TEST(SHAPE_TEST, BESSELT_SET_VARIABLES)
{
    spida::UniformGridRVT grid(4096, -240e-15, 240e-15, 1.1e14, 1.4e16);

    double A = 1.0e8;
    double width = 20.0e-15;
    double apod = 1.0;
    spida::BesselT shape(grid, A, width, apod);

    EXPECT_DOUBLE_EQ(A, shape.amplitude());
    EXPECT_DOUBLE_EQ(width, shape.width());
    EXPECT_DOUBLE_EQ(apod, shape.apodization());
}

TEST(SHAPE_TEST,GAUSST_COMPUTE)
{
    int N = 4096;
    double tmin = -240e-15;
    double tmax = 240e-15;
    double freq_min = 1.10803e14;
    double freq_max = 1.448963e16;
    spida::UniformGridRVT grid(N, tmin, tmax, freq_min, freq_max);

    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;

    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(car_freq);
	auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));
	EXPECT_NEAR(223606797.74997896,*itmax,1e-6);
    EXPECT_NEAR(-200519096.49044815,*itmin,1e-6);
	EXPECT_NEAR(90533935.926272079,y[(N-N/16)/2],1e-6);
	EXPECT_NEAR(90533935.926272079,y[(N+N/16)/2],1e-6);
}

TEST(SHAPE_TEST,SECHT_COMPUTE)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;

    spida::UniformGridRVT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::SechT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(car_freq);
	auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_NEAR(223606797.74997896,*itmax,1e-6);
    EXPECT_NEAR(-211808302.73070601,*itmin,1e-6);
	EXPECT_NEAR(122726544.58500151,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(122726544.58500151,y[(nt+nt/16)/2],1e-6);
}

TEST(SHAPE_TEST,AIRYT_COMPUTE)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -320e-15;
    double maxT = 120e-15;
    double apod = 0.2;

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::AiryT shape(grid,std::sqrt(I0),tp,apod);
    shape.setFastPhase(omega0);
	auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_NEAR(111197223.83176377,*itmax,1e-6);
    EXPECT_NEAR(-115045662.1447591,*itmin,1e-6);
	EXPECT_NEAR(8518347.9235188272,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(-18030265.656660724,y[(nt+nt/16)/2],1e-6);
}

TEST(SHAPE_TEST,BESSELT_COMPUTE)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double omega0 = 4.7091e14;
    double minT = -240e-15;
    double maxT = 240e-15;
    double apod = 0.5;

    spida::UniformGridRVT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::BesselT shape(grid,std::sqrt(I0),tp,apod);
    shape.setFastPhase(omega0);
	auto y = shape.shapeRV();

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_DOUBLE_EQ(223606797.74997896,*itmax);
    EXPECT_DOUBLE_EQ(-185314259.84787187,*itmin);
	EXPECT_DOUBLE_EQ(46643813.615054831,y[(nt-nt/16)/2]);
	EXPECT_DOUBLE_EQ(46643813.615054831,y[(nt+nt/16)/2]);
}