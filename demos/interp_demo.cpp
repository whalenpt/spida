/*------------------------------------------------------------------------------
 *
 *    Author: Patrick Whalen
 *    Email: whalenpt@gmail.com
 *    Status: Completed
 *    Date: 08/17/21
 *    Description: Examples of using interpolation
 *
------------------------------------------------------------------------------*/

#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <spida/helper/interp.h>

static void writeXY(const std::filesystem::path& dir,
                    const std::string& name,
                    const std::vector<double>& x,
                    const std::vector<double>& y)
{
    std::filesystem::create_directories(dir);
    nlohmann::json j{{"type", "xy"}, {"x", x}, {"y", y}};
    std::ofstream(dir / (name + ".json")) << j.dump(2);
}

int main()
{
    std::filesystem::path outdir("outfolder");

    std::vector<double> x(11);
    for(int i = 0; i <= 10; i++)
        x[i] = static_cast<double>(i);
    std::vector<double> y(11);
    for(int i = 0; i <= 10; i++)
        y[i] = sin(x[i]);

    spida::LinearInterp interp(x,y);
    int N = 40;
    std::vector<double> xinterp(N);
    for(int i = 0; i < N; i++)
        xinterp[i] = 10.0*i/(N-1);
    std::vector<double> yinterp = interp.eval(xinterp);

    writeXY(outdir, "data", x, y);
    writeXY(outdir, "interp_data", xinterp, yinterp);

    spida::SplineInterp sinterp(x,y);
    std::vector<double> ysinterp = sinterp.eval(xinterp);
    writeXY(outdir, "spline_interp_data", xinterp, ysinterp);

    // tridisolve
    int n = 20;
    std::vector<double> a(n-1), b(n), c(n-1), d(n);
    for(int i = 0; i < n-1; i++){
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;
        d[i] = i+1;
    }
    b[n-1] = 2.0;
    d[n-1] = n;
    std::vector<double> r;
    spida::tridisolve(a,b,c,d,r);

    std::vector<double> idx(r.size());
    for(std::size_t i = 0; i < r.size(); i++) idx[i] = static_cast<double>(i);
    writeXY(outdir, "tridisolve", idx, r);

    return 0;
}
