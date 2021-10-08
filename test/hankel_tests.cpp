

#include <gtest/gtest.h>
#include <spida/helper/constants.h>
#include <spida/grid/besselR.h>
#include <spida/transform/hankelR.h>
#include <pwutils/report/dataio.hpp>
#include <pwutils/pwmath.hpp>
#include <random>


// Hankel transform of a gauss
// H0{exp(-ar^2)}=(1/2a)*exp(-kr^2/(4a))
// H0{exp(-(r/w0)^2)}=(w0^2/2)*exp(-w0^2*kr^2/4)
TEST(HANKEL_TRANSFORM_TEST,GAUSS)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);

    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    double a = 5.0;
    for(auto i = 0; i < N; i++)
        in[i] = exp(-a*pow(r[i],2));
    for(auto i = 0; i < N; i++){
        double beta = 1.0/(2.0*a);
        exact[i] = beta*exp(-pow(kr[i],2)/(4.0*a));
    }
    transform.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(out,exact),1e-6);
}

// H0{exp(-a*r)/r} = 1.0/(sqrt(a^2+kr^2))
TEST(HANKEL_TRANSFORM_TEST,EXP_OVER_R)
{
    int N = 64;
    double rmax = 6.0;
    double a = 1.0;

    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);
    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    for(auto i = 0; i < N; i++)
        in[i] = exp(-a*r[i])/r[i];
    for(auto i = 0; i < N; i++){
        exact[i] = 1.0/sqrt(pow(a,2)+pow(kr[i],2));
    }
    transform.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(exact,out),0.2);
}

// H0 {sinc(a*r)} = 1.0/(a^2*sqrt(1.0-(kr/a)^2)) if kr < a
//                             0.0               if kr >=a                                        
// sinc(a*r) = sin(a*r)/(a*r)
TEST(HANKEL_TRANSFORM_TEST,SINC_TEST)
{
    int N = 256;
    double rmax = 32.0;
    double a = 5.0;

    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR transform(grid);
    const std::vector<double>& r = grid.getR();
    const std::vector<double>& kr = grid.getSR();

    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> exact(N);

    for(auto i = 0; i < N; i++)
        in[i] = r[i] < 1e-4 ? 1.0-pow(a*r[i],2) : sin(a*r[i])/(a*r[i]);
    for(auto i = 0; i < N; i++)
        exact[i] = kr[i] < a ? 1.0/(pow(a,2)*sqrt(1.0-pow(kr[i]/a,2))) : 0.0;

    transform.R_To_SR(in,out);
    EXPECT_LT(pw::relative_error(exact,out),0.3);
}

// Test that forward Hankel followed by inverse Hankel yields the identity
TEST(HANKEL_TRANSFORM_TEST,INVERSES)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR tr(grid);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,1.0);
    std::vector<double> in(N);
    std::vector<double> out(N);
    std::vector<double> expect(N);
    for(unsigned int i = 0; i < N; i++)
        in[i] = distribution(generator);
    tr.R_To_SR(in,out);
    tr.SR_To_R(out,expect);
    EXPECT_LT(pw::relative_error(in,expect),1e-6);
}

// Test Hankel matrix for orthogonality
TEST(HANKEL_TRANSFORM_TEST,ORTHOGONALITY)
{
    int N = 25;
    double rmax = 2.0;
    spida::BesselRootGridR grid(N,rmax);
    spida::HankelTransformR tr(grid);

    const std::vector<double>& Ymk = tr.getYmk();
    double zero_sum = 0.0;
    for(auto i = 0; i < N; i++){
        for(auto m = 0; m < N; m++){
            double sum = 0.0;
            for(auto k = 0; k < N; k++){
                sum += Ymk[i*N+k]*Ymk[k*N+m];
            }
            if(i == m)
                EXPECT_NEAR(sum,1,1e-6);
            else
                zero_sum += sum;
        }
    }
    EXPECT_NEAR(zero_sum,0,1e-6);
}








