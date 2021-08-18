
#include <gtest/gtest.h>
#include <spida/grid/besselR.h>

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
}



