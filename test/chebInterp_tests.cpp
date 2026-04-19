#include <gtest/gtest.h>
#include <spida/chebInterpX.h>
#include <spida/helper/constants.h>
#include <pwutils/pwmath.hpp>
#include <cmath>
#include <vector>


// xin must cover full [minx,maxx] so the spline can evaluate at interior Chebyshev nodes.
// xout must stay within the Chebyshev root node range (roots don't reach the endpoints).
TEST(CHEB_INTERP_TEST,DERIVATIVE_EXP)
{
    int N_in = 200;
    int N_out = 100;
    double minx = -1.0;
    double maxx = 5.0;

    double dx_in = (maxx - minx) / (N_in - 1);
    std::vector<double> xin(N_in);
    std::vector<double> yin(N_in);
    for(int i = 0; i < N_in; i++){
        xin[i] = minx + i*dx_in;
        yin[i] = exp(xin[i]);
    }

    double margin = 0.02*(maxx - minx);
    double dx_out = (maxx - minx - 2*margin) / (N_out - 1);
    std::vector<double> xout(N_out);
    std::vector<double> expect(N_out);
    std::vector<double> dyout(N_out);
    for(int i = 0; i < N_out; i++){
        xout[i] = minx + margin + i*dx_out;
        expect[i] = exp(xout[i]);
    }

    spida::ChebInterpX interp(64,minx,maxx);
    interp.dXInterp(xin,yin,xout,dyout);
    EXPECT_LT(pw::relative_error(expect,dyout),1e-3);
}

TEST(CHEB_INTERP_TEST,DERIVATIVE_SIN)
{
    using spida::PI;
    int N_in = 200;
    int N_out = 100;
    double minx = 0.0;
    double maxx = PI;

    double dx_in = (maxx - minx) / (N_in - 1);
    std::vector<double> xin(N_in);
    std::vector<double> yin(N_in);
    for(int i = 0; i < N_in; i++){
        xin[i] = minx + i*dx_in;
        yin[i] = sin(xin[i]);
    }

    double margin = 0.02*(maxx - minx);
    double dx_out = (maxx - minx - 2*margin) / (N_out - 1);
    std::vector<double> xout(N_out);
    std::vector<double> expect(N_out);
    std::vector<double> dyout(N_out);
    for(int i = 0; i < N_out; i++){
        xout[i] = minx + margin + i*dx_out;
        expect[i] = cos(xout[i]);
    }

    spida::ChebInterpX interp(64,minx,maxx);
    interp.dXInterp(xin,yin,xout,dyout);
    EXPECT_LT(pw::relative_error(expect,dyout),1e-3);
}

TEST(CHEB_INTERP_TEST,POINT_DERIVATIVE_EXP)
{
    int N_in = 200;
    double minx = -1.0;
    double maxx = 5.0;
    double dx = (maxx - minx) / (N_in - 1);

    std::vector<double> xin(N_in);
    std::vector<double> yin(N_in);
    for(int i = 0; i < N_in; i++){
        xin[i] = minx + i*dx;
        yin[i] = exp(xin[i]);
    }

    spida::ChebInterpX interp(64,minx,maxx);
    double xout = 2.0;
    double dyout = interp.dXInterp(xin,yin,xout);
    EXPECT_NEAR(dyout,exp(xout),1e-3);
}
