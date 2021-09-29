/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Townsend Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/27/21
 *    Description: Implementation of kdv PDE
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/SpidaX.h>
#include <spida/grid/uniformX.h>
#include <spida/helper/constants.h>
#include <spida/solver/solverETDAS.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using namespace spida;
class kDV
{
    public: 
        kDV(const UniformGridX& grid) : 
            L(grid.getNsx()),
            m_grid(grid), 
            m_spida(grid), 
            m_uphys(grid.getNx()),
            m_uxphys(grid.getNx()),
            m_uxsp(grid.getNsx())
            {
                const std::vector<double>& sx = grid.getSX();
                for(auto i = 0; i < sx.size(); i++)
                    L[i] = ii*pow(sx[i],3);
            }
        std::vector<dcmplx> L;
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> NL = [this](\
                const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
            m_spida.SX_To_X(in,m_uphys);
            m_spida.dSX(in,m_uxsp);
            m_spida.SX_To_X(m_uxsp,m_uxphys);
            for(auto i = 0; i < m_grid.getNx(); i++)
                m_uphys[i] = -6.0*m_uphys[i]*m_uxphys[i];
            m_spida.X_To_SX(m_uphys,out);
        };
        SpidaX& spida() {return m_spida;}

    private:
        UniformGridX m_grid;
        SpidaX m_spida;
        std::vector<dcmplx> m_uphys;
        std::vector<dcmplx> m_uxphys;
        std::vector<dcmplx> m_uxsp;
};

int main()
{
    unsigned nx = 512;
    double minx = -150.0;
    double maxx = 150.0;
    std::string outdir("kdv_files");

    UniformGridX grid(nx,minx,maxx);
    kDV model(grid);
    ETD34 solver(model.L,model.NL);
    solver.setEpsRel(1e-4);

    std::vector<dcmplx> u0(nx,0.0);
    std::vector<double> A0{0.6,0.5,0.4,0.3,0.2};
    std::vector<double> x0{-120.0,-90.0,-60.0,-30.0,0.0};
    const std::vector<double>& x  = grid.getX();
    for(auto i = 0; i < nx; i++){
        for(auto k = 0; k < A0.size(); k++){
            u0[i] += 0.5*pow(A0[k],2)/pow(cosh(A0[k]*(x[i]-x0[k])/2.0),2);
        }
    }
    std::vector<dcmplx> usp(grid.getNsx());
    model.spida().X_To_SX(u0,usp);

    double t0 = 0.0;
    double tf = 600.0;
    double h = 0.1;
    double h_next = 0.1;
    double t = t0;

    std::vector<dcmplx> uphys(nx);
    std::copy(std::cbegin(u0),std::cend(u0),std::begin(uphys));
    dat::ReportData1D<double,dcmplx> report("X",x,uphys);
    report.setDirPath(outdir);
    report.setItem("t",t0);
    std::cout << "First physical space report file location: " << report.path(0) << std::endl;

    std::vector<double> shifted_kx = grid.freqshift(grid.getSX());
    std::vector<dcmplx> shifted_usp = grid.freqshift(usp);
    dat::ReportData1D<double,dcmplx> reportS("SX",shifted_kx,shifted_usp);
    reportS.setDirPath(outdir);
    std::cout << "First spectral space report file location: " << reportS.path(0) << std::endl;

    std::ofstream os;
    os.open(report.path(0));
    os << std::scientific << std::setprecision(8);
    report.report(os);
    os.close();

    os.open(reportS.path(0));
    reportS.report(os);
    os.close();

    unsigned step_count = 0;
    unsigned report_count = 1;
    while(t < tf){
        if(!solver.step(usp,h,h_next)){
            std::cout << "Step failed." << std::endl;
            return 1;
        }
        t += h;
        h = h_next;
        model.spida().SX_To_X(usp,uphys);
        if(step_count % 4 == 0){
            os.open(report.path(report_count));
            report.report(os);
            report.setItem("t",t);
            os.close();
            grid.freqshift(usp,shifted_usp);
            os.open(reportS.path(report_count));
            reportS.report(os);
            os.close();
            report_count++;
        }
        step_count++;
    }
    os.open(report.path(report_count));
    report.report(os);
    report.setItem("t",t);
    os.close();
    grid.freqshift(usp,shifted_usp);
    os.open(reportS.path(report_count));
    reportS.report(os);
    os.close();

    return 0;
}





