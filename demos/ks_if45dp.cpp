/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Complete
 *    Date: 10/24/21
 *    Description: Implementation of KS PDE with Dormand-Prince 
 *                 integrating factor method IF45DP
 *                   
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/RVX.h>
#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/IFAS.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using namespace spida;

// KS model for real-valued physical space fields (spectral space is complex)
class KS_RV
{
    public: 
        KS_RV(const UniformGridRVX& grid) : 
            m_grid(grid), 
            m_spi(grid), 
            m_uphys(grid.getNx()),
            m_uxphys(grid.getNx()),
            m_uxsp(grid.getNsx()),
            m_L(grid.getNsx())
            {
                const std::vector<double>& sx = grid.getSX();
                for(auto i = 0; i < sx.size(); i++)
                    m_L[i] = pow(sx[i],2)*(1.0-pow(sx[i],2));
                m_NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    m_spi.SX_To_X(in,m_uphys);
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
class PropagatorKS : public PropagatorCV
{
    public:
        PropagatorKS(const std::filesystem::path& path,KS_RV& md) : 
            PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().getGridX().getNsx(),0.0),
            m_uphys(md.spida().getGridX().getNx(),0.0) 
         {
             // initialize propagator m_usp
             const std::vector<double>& x  = m_spi.getX();
             for(auto i = 0; i < x.size(); i++)
                 m_uphys[i] = cos(x[i]/16.0)*(1.0+sin(x[i]/16.0));
             // Need to initialize the propagator which is the spectral space representation of m_uphys
             m_spi.X_To_SX(m_uphys,m_usp);
             initReport();
         }
        ~PropagatorKS() {}
        // updateFields is a pure virtual function of PropagatorCV and must be implemented 
        // This function is called before each Solver report (allows for updating of real space fields)
        void updateFields(double t) { m_spi.SX_To_X(m_usp,m_uphys);}
        std::vector<dcmplx>& propagator() {return m_usp;}
    private:
        // initReport is a helper function that feeds PropagatorCV information on what to report out to files
        void initReport() {
            // add report for real space KS field
            const std::vector<double>& x  = m_spi.getGridX().getX();
            PropagatorCV::addReport(std::make_unique<dat::ReportData1D<double,double>>("X",x,m_uphys));
            // add report for spectral space KS field (the propagator)
            const std::vector<double>& sx  = m_spi.getGridX().getSX();
            PropagatorCV::addReport(std::make_unique<dat::ReportComplexData1D<double,double>>("SX",sx,m_usp));
        }
        SpidaRVX& m_spi;
        std::vector<dcmplx> m_usp;
        std::vector<double> m_uphys;
};

int main()
{
    unsigned N = 8192;
    double a = 0.0;
    double b = 32.0*spida::PI;

    UniformGridRVX grid(N,a,b);
    KS_RV model(grid);

    std::filesystem::path dirpath("ksif45dp_propagator_files");
    PropagatorKS propagator(dirpath,model);
    propagator.setStepsPerOutput(15);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(200);
    // End simulation once 100 files have been output
    propagator.setMaxReports(100);

    IF45DP solver(model.L(),model.NL());
    solver.setEpsRel(1e-3);
    solver.setLogProgress(true);
    solver.setLogFrequency(200);
    solver.evolve(propagator,0.0,50.0,0.5);

    return 0;
}





