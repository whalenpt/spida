/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 10/24/21
 *    Description: Implementation of KS PDE with an integrating factor method
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

using spida::dcmplx;
using spida::ii;

// KS model for real-valued physical space fields (spectral space is complex)
class KS_RV
{
    public: 
        explicit KS_RV(const spida::UniformGridRVX& grid) : 
            m_grid(grid), 
            m_spi(grid), 
            m_uphys(grid.getNx()),
            m_uxphys(grid.getNx()),
            m_uxsp(grid.getNsx()),
            m_L(grid.getNsx())
            {
                const auto& sx = grid.getSX();
                for(size_t i = 0; i < sx.size(); i++)
                    m_L[i] = pow(sx[i],2)*(1.0-pow(sx[i],2));
                m_NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    m_spi.SX_To_X(in,m_uphys);
                    m_spi.dSX(in,m_uxsp);
                    m_spi.SX_To_X(m_uxsp,m_uxphys);
                    for(unsigned i = 0; i < m_grid.getNx(); i++)
                        m_uphys[i] = -m_uphys[i]*m_uxphys[i];
                    m_spi.X_To_SX(m_uphys,out);
                };
            }
        std::vector<dcmplx>& L() {return m_L;}
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)>& NL() {return m_NL;}
        spida::SpidaRVX& spida() {return m_spi;}

    private:
        spida::UniformGridRVX m_grid;
        spida::SpidaRVX m_spi;
        std::vector<double> m_uphys;
        std::vector<double> m_uxphys;
        std::vector<dcmplx> m_uxsp;
        std::vector<dcmplx> m_L;
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> m_NL;
};

// Helper class for reporting files based on data generated from the Solver used
class PropagatorKS : public spida::PropagatorCV
{
    public:
        PropagatorKS(const std::filesystem::path& path,KS_RV& md) : 
            PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().getGridX().getNsx(),0.0),
            m_uphys(md.spida().getGridX().getNx(),0.0) 
         {
             // initialize propagator m_usp
             const auto& x  = m_spi.getX();
             for(size_t i = 0; i < x.size(); i++)
                 m_uphys[i] = cos(x[i]/16.0)*(1.0+sin(x[i]/16.0));
             // Need to initialize the propagator which is the spectral space representation of m_uphys
             m_spi.X_To_SX(m_uphys,m_usp);
             initReport();
         }
        ~PropagatorKS() override = default;
        // updateFields is a pure virtual function of PropagatorCV and must be implemented 
        // This function is called before each Solver report (allows for updating of real space fields)
        void updateFields(double t) override { m_spi.SX_To_X(m_usp,m_uphys);}
        std::vector<dcmplx>& propagator() override {return m_usp;}
    private:
        // initReport is a helper function that feeds PropagatorCV information on what to report out to files
        void initReport() {
            // add report for real space KS field
            const auto& x  = m_spi.getGridX().getX();
            PropagatorCV::addReport(std::make_unique<dat::ReportData1D<double,double>>("X",x,m_uphys));
            // add report for spectral space KS field (the propagator)
            const auto& sx  = m_spi.getGridX().getSX();
            PropagatorCV::addReport(std::make_unique<dat::ReportComplexData1D<double,double>>("SX",sx,m_usp));
        }
        spida::SpidaRVX& m_spi;
        std::vector<dcmplx> m_usp;
        std::vector<double> m_uphys;
};

int main()
{
    unsigned N = 8192;
    double a = 0.0;
    double b = 32.0*spida::PI;

    spida::UniformGridRVX grid(N,a,b);
    KS_RV model(grid);

    std::filesystem::path dirpath("ksif34_propagator_files");
    PropagatorKS propagator(dirpath,model);
    propagator.setStepsPerOutput(5);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(200);

    spida::IF34 solver(model.L(),model.NL());
    solver.setEpsRel(1e-4);
    solver.setLogProgress(true);
    solver.setLogFrequency(200);
    solver.evolve(propagator,0.0,50.0,0.5);

    return 0;
}





