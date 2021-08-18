
#include <gtest/gtest.h>
#include <spida/grid/besselR.h>
#include <boost/math/special_functions/bessel.hpp>

TEST(GRID_TEST,BESSELZERO_R)
{
    int N = 5;
    double accuracy = 1e-4;
    spida::BesselRootGridR grid = spida::BesselRootGridR(N,1.0);
    const std::vector<double>& roots = grid.getBesselRoots();
	EXPECT_NEAR(roots[0],2.4048,accuracy);
	EXPECT_NEAR(roots[1],5.5201,accuracy);
	EXPECT_NEAR(roots[2],8.6537,accuracy);
	EXPECT_NEAR(roots[3],11.7915,accuracy);
	EXPECT_NEAR(roots[4],14.9309,accuracy);

    std::vector<double> J1;
    boost::math::cyl_bessel_j_zero(1.0,1,N,std::back_inserter(J1));
	EXPECT_NEAR(J1[0],3.8317,accuracy);
	EXPECT_NEAR(J1[1],7.0156,accuracy);
	EXPECT_NEAR(J1[2],10.1735,accuracy);
	EXPECT_NEAR(J1[3],13.3237,accuracy);
	EXPECT_NEAR(J1[4],16.4706,accuracy);

}



