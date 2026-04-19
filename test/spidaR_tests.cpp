#include <gtest/gtest.h>
#include <spida/R.h>
#include <spida/helper/constants.h>
#include <spida/grid/besselR.h>
#include <pwutils/pwmath.hpp>
#include <random>


TEST(SPIDA_R_TEST,INVERSES)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::SpidaR spidaR(grid);

    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);
    for(int i = 0; i < N; i++)
        in[i] = distribution(generator);
    spidaR.R_To_SR(in,out);
    spidaR.SR_To_R(out,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

TEST(SPIDA_R_TEST,INVERSES_COMPLEX)
{
    using spida::dcmplx;
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::SpidaR spidaR(grid);

    std::default_random_engine generator;
    std::normal_distribution distribution{1.0,1.0};
    std::vector<dcmplx> in(N);
    std::vector<dcmplx> out(N);
    std::vector<dcmplx> expect(N);
    for(int i = 0; i < N; i++)
        in[i] = dcmplx(distribution(generator),distribution(generator));
    spidaR.R_To_SR(in,out);
    spidaR.SR_To_R(out,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

// H0{exp(-ar^2)}=(1/2a)*exp(-kr^2/(4a))
TEST(SPIDA_R_TEST,GAUSS)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::SpidaR spidaR(grid);

    const auto& r = grid.getR();
    const auto& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    double a = 5.0;
    for(int i = 0; i < N; i++)
        in[i] = exp(-a*pow(r[i],2));
    for(int i = 0; i < N; i++){
        double beta = 1.0/(2.0*a);
        exact[i] = beta*exp(-pow(kr[i],2)/(4.0*a));
    }
    spidaR.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);
}
