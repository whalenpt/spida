/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Townsend Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/27/21
 *    Description: Implementation of kdv PDE using complex physical space values
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/uniformCVX.h>
#include <spida/SpidaCVX.h>
#include <spida/helper/constants.h>
#include <spida/rkstiff/ETDAS.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using namespace spida;
class kDV
{
    public: 
        kDV(const UniformGridCVX& grid) : 
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
        SpidaCVX& spida() {return m_spida;}

    private:
        UniformGridCVX m_grid;
        SpidaCVX m_spida;
        std::vector<dcmplx> m_uphys;
        std::vector<dcmplx> m_uxphys;
        std::vector<dcmplx> m_uxsp;
};

int main()
{
    unsigned nx = 512;
    double minx = -150.0;
    double maxx = 150.0;
    std::string outdir("kdv_files_CV");

    UniformGridCVX grid(nx,minx,maxx);
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
    dat::ReportData1D<double,dcmplx> report("X_0",x,uphys);
    report.setDirPath(outdir);
    report.setItem("t",t0);
    std::cout << "First physical space report file location: " << report.path() << std::endl;

    std::vector<double> shifted_kx = grid.freqshift(grid.getSX());
    std::vector<dcmplx> shifted_usp = grid.freqshift(usp);
    dat::ReportData1D<double,dcmplx> reportS("SX_0",shifted_kx,shifted_usp);
    reportS.setDirPath(outdir);
    std::cout << "First spectral space report file location: " << reportS.path() << std::endl;

    std::ofstream os;
    os << std::scientific << std::setprecision(8);
    os << report;
    os << reportS;

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
        if(step_count % 16 == 0){
            report.setName("X_" + std::to_string(report_count));
            reportS.setName("SX_" + std::to_string(report_count));
            report.setItem("t",t);
            reportS.setItem("t",t);
            grid.freqshift(usp,shifted_usp);
            os << report;
            os << reportS;
            report_count++;
        }
        step_count++;
    }
    report.setName("X_" + std::to_string(report_count));
    reportS.setName("SX_" + std::to_string(report_count));
    report.setItem("t",t);
    reportS.setItem("t",t);
    os << report;
    os << reportS;
    os.close();
    return 0;
}





