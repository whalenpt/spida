/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Townsend Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/27/21
 *    Description: Implementation of Burgers PDE with a Propagator(
 *    automated file reporting)
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/uniformRVX.h>
#include <spida/SpidaRVX.h>
#include <spida/helper/constants.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/propagator/propagator.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using namespace spida;

// Burgers model for real-valued physical space fields (spectral space is complex)
class Burgers
{
    public: 
        Burgers(const UniformGridRVX& grid) : 
            m_grid(grid), 
            m_spi(grid), 
            m_uphys(grid.getNx()),
            m_uxphys(grid.getNx()),
            m_uxsp(grid.getNsx()),
            m_L(grid.getNsx())
            {
                double mu = 0.0005;
                const std::vector<double>& sx = grid.getSX();
                for(auto i = 0; i < sx.size(); i++)
                    m_L[i] = -mu*pow(sx[i],2);
                m_NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    m_spi.SX_To_X(in,m_uphys);
                    // dSX takes derivative in spectral space, i.e. multiply by i*kx
                    m_spi.dSX(in,m_uxsp);
                    m_spi.SX_To_X(m_uxsp,m_uxphys);
                    for(auto i = 0; i < m_grid.getNx(); i++)
                        m_uphys[i] = -m_uphys[i]*m_uxphys[i];
                    m_spi.X_To_SX(m_uphys,out);
                };
            }
        std::vector<dcmplx>& L() {return m_L;}
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)>& NL() {return m_NL;}
        SpidaRVX& spida() {return m_spi;}

    private:
        UniformGridRVX m_grid;
        SpidaRVX m_spi;
        std::vector<double> m_uphys;
        std::vector<double> m_uxphys;
        std::vector<dcmplx> m_uxsp;
        std::vector<dcmplx> m_L;
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> m_NL;
};

// Helper class for reporting files based on data generated from the Solver used
class BurgProp : public PropagatorCV
{
    public:
        BurgProp(const std::filesystem::path& path,Burgers& md) : 
            PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().getGridX().getNsx(),0.0),
            m_uphys(md.spida().getGridX().getNx(),0.0) 
         {
             // initialize propagator m_usp
             const std::vector<double>& x  = m_spi.getX();
             for(auto i = 0; i < x.size(); i++)
                 m_uphys[i] = exp(-10.0*pow(sin(x[i]/2.0),2));
             // Need to initialize the propagator which is the spectral space representation of m_uphys
             m_spi.X_To_SX(m_uphys,m_usp);
             // initialize reporting arrays
             initReport();
         }
        ~BurgProp() {}
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
    // construct grids
    unsigned N = 8192;
    double a = -spida::PI;
    double b = spida::PI;
    UniformGridRVX grid(N,a,b);

    // construct linear operator and nonlinear function
    Burgers model(grid);

    // file reporting aide
    std::filesystem::path dirpath("burgers_files");
    BurgProp propagator(dirpath,model);
    propagator.setStepsPerOutput(5);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(200);

    // rkstiff solver
    ETD35 solver(model.L(),model.NL());
    solver.setEpsRel(1e-8);
    solver.setLogProgress(true);
    solver.setLogFrequency(200);
    solver.setModeCutoff(0.001);
    solver.evolve(propagator,0.0,1.0,0.5);

    return 0;
}





