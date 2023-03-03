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

#include <spida/R.h>
#include <spida/grid/besselR.h>
#include <spida/helper/constants.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/propagator/propagator.h>
#include <pwutils/report/dat.hpp>
#include <complex> // std::norm -> std::norm(dcmplx(3,4)) = 25

//------------------------------------------------------------------------------

using namespace spida;

// NLS model for complex-valued physical space fields and spectral fields
class NLSR
{
    public: 
        explicit NLSR(const BesselRootGridR& grid) : 
            m_spi(grid), 
            m_uphys(grid.getNr()),
            m_L(grid.getNsr())
            {
                const std::vector<double>& kr = grid.getSR();
                // ii is spida::ii which is imaginary an number
                for(size_t i = 0; i < kr.size(); i++)
                   m_L[i] = -ii*pow(kr[i],2);
                 m_NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    double gamma = 2.0;
                    m_spi.SR_To_R(in,m_uphys);
                    for(auto& val : m_uphys)
                        val = ii*gamma*val*std::norm(val);
                    m_spi.R_To_SR(m_uphys,out);
                };
            }
        std::vector<dcmplx>& L() {return m_L;}
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)>& NL() {return m_NL;}
        SpidaR& spida() {return m_spi;}

    private:
        SpidaR m_spi;
        std::vector<dcmplx> m_uphys;
        std::vector<dcmplx> m_L;
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> m_NL;
};

// Helper class for reporting files based on data generated from the Solver used
class PropagatorNLSR : public PropagatorCV
{
    public:
        PropagatorNLSR(const std::filesystem::path& path,NLSR& md) : 
            PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().getGridR().getNsr(),0.0),
            m_uphys(md.spida().getGridR().getNr(),0.0),
            m_mirror_r(2*md.spida().getGridR().getNr()),
            m_mirror_uphys(2*m_uphys.size()),
            m_mirror_kr(2*md.spida().getGridR().getNsr()),
            m_mirror_usp(2*m_usp.size())
         {
             // initialize propagator m_usp
             const std::vector<double>& r  = m_spi.getR();
             double A0 = 2.0; // amplitude
             double a = 1.0; // width of gaussian
             for(size_t i = 0; i < r.size(); i++)
                 m_uphys[i] = A0*exp(-a*pow(r[i],2));
             // Need to initialize the propagator which is the spectral space representation of m_uphys
             m_spi.R_To_SR(m_uphys,m_usp);
             // mirrored radial components for better graphs
             m_mirror_r = m_spi.getGridR().mirrorGrid(m_spi.getR(),true);
             // mirrored spectral components for better graphs
             m_mirror_kr = m_spi.getGridR().mirrorGrid(m_spi.getSR(),true);
             // mirrored physical space field
             m_spi.getGridR().mirrorGrid(m_uphys,m_mirror_uphys);
             // mirrored spectral space field
             m_spi.getGridR().mirrorGrid(m_usp,m_mirror_usp);
             initReport();
         }
        ~PropagatorNLSR() override = default;
        // updateFields is a pure virtual function of PropagatorCV and must be implemented 
        // This function is called before each Solver report (allows for updating of real space fields)
        void updateFields(double t) override { 
            m_spi.SR_To_R(m_usp,m_uphys);
            // output mirrored grids with negative radial components (for better graphs, not necessary)
            m_spi.getGridR().mirrorGrid(m_uphys,m_mirror_uphys);
            m_spi.getGridR().mirrorGrid(m_usp,m_mirror_usp);
        }
        std::vector<dcmplx>& propagator() override {return m_usp;}
    private:
        // initReport is a helper function that feeds PropagatorCV information on what to report out to files
        void initReport() {
            // add report for complex physical space NLS field
            auto report = std::make_unique<dat::ReportComplexData1D<double,double>>("R",m_mirror_r,m_mirror_uphys);
            PropagatorCV::addReport(std::move(report));

            // add report for power of physical space NLS field 
            auto reportpow = std::make_unique<dat::ReportComplexData1D<double,double>>("SQ_R",m_mirror_r,m_mirror_uphys);
            reportpow->setPower(true);
            PropagatorCV::addReport(std::move(reportpow));

            // add report for spectral space NLS field (the propagator)
            auto reportsp = std::make_unique<dat::ReportComplexData1D<double,double>>("SR",m_mirror_kr,m_mirror_usp);
            PropagatorCV::addReport(std::move(reportsp));

            // add report for power of spectral space NLS field (the propagator)
            auto reportsppow = std::make_unique<dat::ReportComplexData1D<double,double>>("SQ_SR",m_mirror_kr,m_mirror_usp);
            reportsppow->setPower(true);
            PropagatorCV::addReport(std::move(reportsppow));
        }
        SpidaR& m_spi;
        std::vector<dcmplx> m_usp;
        std::vector<dcmplx> m_uphys;
        std::vector<double> m_mirror_r;
        std::vector<dcmplx> m_mirror_uphys;
        std::vector<double> m_mirror_kr;
        std::vector<dcmplx> m_mirror_usp;
};

int main()
{
    unsigned N = 100;
    double rmax = 5.0;
    BesselRootGridR grid(N,rmax);
    NLSR model(grid);

    std::filesystem::path dirpath("nlsR_propagator_files");
    PropagatorNLSR propagator(dirpath,model);
    propagator.setStepsPerOutput(10);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(10);

    ETD35 solver(model.L(),model.NL());
    solver.setEpsRel(1e-5);
    solver.setLogProgress(true);
    solver.setLogFrequency(10);
    solver.evolve(propagator,0.0,0.8,0.01);
    return 0;
}