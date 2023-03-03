/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Townsend Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/27/21
 *    Description: Implementation of kdv PDE with a Propagator class (automated file reporting)
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/RVX.h>
#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/ETDAS.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using namespace spida;

// kDV model for real-valued physical space fields (spectral space is still complex)
class kDV_RV
{
    public: 
        kDV_RV(const UniformGridRVX& grid) : 
            m_grid(grid), 
            m_spida(grid), 
            m_uphys(grid.getNx()),
            m_uxphys(grid.getNx()),
            m_uxsp(grid.getNsx()),
            m_L(grid.getNsx())
            {
                const std::vector<double>& sx = grid.getSX();
                for(auto i = 0; i < sx.size(); i++)
                    m_L[i] = ii*pow(sx[i],3);
                m_NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    m_spida.SX_To_X(in,m_uphys);
                    m_spida.dSX(in,m_uxsp);
                    m_spida.SX_To_X(m_uxsp,m_uxphys);
                    for(auto i = 0; i < m_grid.getNx(); i++)
                        m_uphys[i] = -6.0*m_uphys[i]*m_uxphys[i];
                    m_spida.X_To_SX(m_uphys,out);
                };
            }
        std::vector<dcmplx>& L() {return m_L;}
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)>& NL() {return m_NL;}
        SpidaRVX& spida() {return m_spida;}

    private:
        UniformGridRVX m_grid;
        SpidaRVX m_spida;
        std::vector<double> m_uphys;
        std::vector<double> m_uxphys;
        std::vector<dcmplx> m_uxsp;
        std::vector<dcmplx> m_L;
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> m_NL;
};

// Helper class for reporting files based on data generated from the Solver used
class PropagatorKDV : public PropagatorCV
{
    public:
        PropagatorKDV(const std::filesystem::path& path,kDV_RV& md) : PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().getGridX().getNsx(),0.0),
            m_uphys(md.spida().getGridX().getNx(),0.0) 
         {
             // initialize propagator m_usp
             std::vector<double> A0{0.6,0.5,0.4,0.3,0.2};
             std::vector<double> x0{-120.0,-90.0,-60.0,-30.0,0.0};
             const std::vector<double>& x  = m_spi.getX();
             for(auto i = 0; i < x.size(); i++){
                 for(auto k = 0; k < A0.size(); k++){
                     m_uphys[i] += 0.5*pow(A0[k],2)/pow(cosh(A0[k]*(x[i]-x0[k])/2.0),2);
                 }
             }
             // Need to initialize the propagator which is the spectral space representation of m_uphys
             m_spi.X_To_SX(m_uphys,m_usp);
             initReport();
         }
        ~PropagatorKDV() {}
        // updateFields is a pure virtual function of PropagatorCV and must be implemented 
        // This function is called before each Solver report (allows for updating of real space fields)
        void updateFields(double t) { m_spi.SX_To_X(m_usp,m_uphys);}
        std::vector<dcmplx>& propagator() {return m_usp;}
    private:
        // initReport is a helper function that feeds PropagatorCV information on what to report out to files
        void initReport() {
            // add report for real space kDV field
            const std::vector<double>& x  = m_spi.getGridX().getX();
            auto report = std::make_unique<dat::ReportData1D<double,double>>("X",x,m_uphys);
            PropagatorCV::addReport(std::move(report));
            // add report for spectral space kDV field (the propagator)
            const std::vector<double>& sx  = m_spi.getGridX().getSX();
            auto reportsp = std::make_unique<dat::ReportComplexData1D<double,double>>("SX",sx,m_usp);
            PropagatorCV::addReport(std::move(reportsp));
        }
        SpidaRVX& m_spi;
        std::vector<dcmplx> m_usp;
        std::vector<double> m_uphys;
};

int main()
{
    unsigned nx = 512;
    double minx = -150.0;
    double maxx = 150.0;

    UniformGridRVX grid(nx,minx,maxx);
    kDV_RV model(grid);

    std::filesystem::path dirpath("kdv_propagator_files");
    PropagatorKDV propagator(dirpath,model);
    propagator.setStepsPerOutput(16);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(200);

    ETD34 solver(model.L(),model.NL());
    solver.setEpsRel(1e-4);
    solver.setLogProgress(true);
    solver.setLogFrequency(200);
    solver.evolve(propagator,0.0,600.0,0.1);

    return 0;
}





