
#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/uniformT.h>
#include <spida/shape/shapeT.h>
#include <pwutils/report/dataio.hpp>

TEST(SHAPE_TEST,GAUSST)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;
    std::vector<double> y(nt);

    spida::UniformGridT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::GaussT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(car_freq);
    shape.shapeRV(y);


    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));
	EXPECT_NEAR(223519721.79339397,*itmax,1e-6);
    EXPECT_NEAR(-200543934.79216298,*itmin,1e-6);
	EXPECT_NEAR(93205462.4168476,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(87456376.64073436,y[(nt+nt/16)/2],1e-6);
}

TEST(SHAPE_TEST,SECHT)
{
    int nt = 4096;
    double I0 = 5.0e16;
    double tp = 20.0e-15;
    double car_freq = 4.7091e14;
    std::vector<double> y(nt);

    spida::UniformGridT grid(nt,-240e-15,240e-15,1.10803e14,1.448963e16);
    spida::SechT shape(grid,std::sqrt(I0),tp);
    shape.setFastPhase(car_freq);
    shape.shapeRV(y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_NEAR(223520681.50796005,*itmax,1e-6);
    EXPECT_NEAR(-211828145.78948778,*itmin,1e-6);
	EXPECT_NEAR(126048853.90398246,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(118875126.49027915,y[(nt+nt/16)/2],1e-6);
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

    spida::UniformGridT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::AiryT shape(grid,std::sqrt(I0),tp,apod);
    shape.setFastPhase(omega0);
    shape.shapeRV(y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_NEAR(111169075.4861387,*itmax,1e-6);
    EXPECT_NEAR(-115056619.97539623,*itmin,1e-6);
	EXPECT_NEAR(8435720.383434976,y[(nt-nt/16)/2],1e-6);
	EXPECT_NEAR(-17718946.86783132,y[(nt+nt/16)/2],1e-6);
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

    spida::UniformGridT grid(nt,minT,maxT,1.10803e14,1.448963e16);
    spida::BesselT shape(grid,std::sqrt(I0),tp,apod);
    shape.setFastPhase(omega0);
    shape.shapeRV(y);

    auto itmax = std::max_element(std::begin(y),std::end(y));
    auto itmin = std::min_element(std::begin(y),std::end(y));

	EXPECT_DOUBLE_EQ(223518386.26149443,*itmax);
    EXPECT_DOUBLE_EQ(-185330357.5625744,*itmin);
	EXPECT_DOUBLE_EQ(48416628.206692584,y[(nt-nt/16)/2]);
	EXPECT_DOUBLE_EQ(44633862.627727090,y[(nt+nt/16)/2]);
}




    







