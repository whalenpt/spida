/*------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/27/21
 *    Description: Implementation of NLS PDE with a Propagator class (automated file reporting)
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/grid/besselR.h>
#include <spida/grid/uniformCVT.h>
#include <spida/SpidaRCVT.h>
#include <spida/helper/constants.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/propagator/propagator.h>
#include <pwutils/report/dat.hpp>
#include <algorithm>
#include <complex> // std::norm -> std::norm(dcmplx(3,4)) = 25

//------------------------------------------------------------------------------

using namespace spida;
// dzA = -0.5*i*dt^2A -i\laplace_perpA  + i|U|^2U

// NLS model for complex-valued physical space fields and spectral fields
class NLSRT
{
    public: 
        NLSRT(const BesselRootGridR& gridR,const UniformGridCVT& gridT,unsigned threads) : 
            m_spi(gridR,gridT,threads), 
            m_uphys(gridR.getNr()*gridT.getNt()),
            m_L(gridR.getNsr()*gridT.getNst())
            {
                const std::vector<double>& kr = gridR.getSR();
                const std::vector<double>& omega = gridT.getST();
                double sigma = 0.5;
                // ii is spida::ii which is imaginary an number sqrt(-1)
                for(auto i = 0; i < kr.size(); i++)
                    for(auto j = 0; j < omega.size(); j++)
                        m_L[i*omega.size()+j] = -ii*pow(kr[i],2) + ii*sigma*pow(omega[j],2);

                m_NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    double gamma = 2.0;
                    m_spi.SRST_To_RT(in,m_uphys);
                    for(auto i = 0; i < m_uphys.size(); i++)
                        m_uphys[i] = ii*gamma*m_uphys[i]*std::norm(m_uphys[i]);
                    m_spi.RT_To_SRST(m_uphys,out);
                };
            }
        std::vector<dcmplx>& L() {return m_L;}
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)>& NL() {return m_NL;}
        SpidaRCVT& spida() {return m_spi;}

    private:
        SpidaRCVT m_spi;
        std::vector<dcmplx> m_uphys;
        std::vector<dcmplx> m_L;
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> m_NL;
};

// Helper class for reporting files based on data generated from the Solver used
class Propagator : public PropagatorCV
{
    public:
        Propagator(const std::filesystem::path& path,NLSRT& md) : 
            PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().spectralSize(),0.0),
            m_uphys(md.spida().physicalSize(),0.0),
            m_mirror_r(2*md.spida().getGridR().getNr()),
            m_mirror_kr(2*md.spida().getGridR().getNsr()),
            m_shifted_omega(md.spida().getGridT().getNst()),
            m_mirror_uphys(2*m_uphys.size()),
            m_mirror_shift_usp(2*m_usp.size())
         {
             // initialize propagator m_usp
             double A0 = 4.0; // amplitude
             double w0 = 1.0;
             double tp = 0.5; // width of gaussian
             const std::vector<double>& r  = m_spi.getGridR().getR();
             const std::vector<double>& t  = m_spi.getGridT().getT();
             for(auto i = 0; i < r.size(); i++)
                 for(auto j = 0; j < t.size(); j++)
                     m_uphys[i*t.size()+j] = A0*exp(-pow(r[i]/w0,2)-pow(t[j]/tp,2));
             // Need to initialize the propagator which is the spectral space representation of m_uphys
             m_spi.RT_To_SRST(m_uphys,m_usp);

             // mirrored radial components for better graphs
             m_mirror_r = m_spi.getGridR().mirrorGrid(m_spi.getR(),true);

             // mirrored spectral components for better graphs
             m_mirror_kr = m_spi.getGridR().mirrorGrid(m_spi.getSR(),true);

             // mirrored physical space field via R coordinate
             m_spi.mirrorR(m_uphys,m_mirror_uphys);

             // mirrored spectral space field
             // m_spi.getGridR().mirrorGrid(m_usp,m_mirror_usp);
             initReport();
         }
        ~Propagator() {}
        // updateFields is a pure virtual function of PropagatorCV and must be implemented 
        // This function is called before each Solver report (allows for updating of real space fields)
        void updateFields(double t) { 
            //std::cout << "updateFields called" << std::endl;
            m_spi.SRST_To_RT(m_usp,m_uphys);
            m_spi.mirrorR(m_uphys,m_mirror_uphys);
            // output mirrored grids with negative radial components (for better graphs)
            //m_spi.getGridR().mirrorGrid(m_uphys,m_mirror_uphys);
            //m_spi.getGridR().mirrorGrid(m_usp,m_mirror_usp);
        }
        // propagator accessed by solver (what is propagated)
        std::vector<dcmplx>& propagator() {return m_usp;}
    private:
        // initReport is a helper function that feeds PropagatorCV information on what to report out to files
        void initReport() {
            // add report for complex physical space NLS field
            auto report = std::make_unique<dat::ReportComplexData2D<\
                     double,double,double>>("RT",m_mirror_r,m_spi.getT(),m_mirror_uphys);
            report->setItem("xlabel","r");
            report->setItem("ylabel","t");
            report->setItem("zlabel","A");
            report->setStrideX(2);
            report->setStrideY(8);
            PropagatorCV::addReport(std::move(report));

//            auto track_max = std::make_unique<dat::TrackComplexData<double>>("Peak Power",\
//                    pw::TrackType::Max,m_uphys,pw::ComplexOp::Power);
//            PropagatorCV::addReport(std::move(track_max));

            // add report for power of physical space NLS field 
            auto reportpow = std::make_unique<dat::ReportComplexData2D<\
                             double,double,double>>("SQ_RT",m_mirror_r,m_spi.getT(),m_mirror_uphys);
            reportpow->setPower(true);
            reportpow->setItem("xlabel","r");
            reportpow->setItem("ylabel","t");
            reportpow->setItem("zlabel","|A|^2");
            reportpow->setStrideX(2);
            reportpow->setStrideY(8);
            PropagatorCV::addReport(std::move(reportpow));

            // add report for spectral space NLS field (the propagator)
            auto reportsp = std::make_unique<dat::ReportComplexData2D<\
                        double,double,double>>("SR",m_spi.getSR(),m_spi.getST(),m_usp);
            reportsp->setItem("xlabel","kr");
            reportsp->setItem("ylabel","omega");
            reportsp->setItem("zlabel","A");
            reportsp->setStrideX(2);
            reportsp->setStrideY(8);
            PropagatorCV::addReport(std::move(reportsp));

//            // add report for power of spectral space NLS field (the propagator)
//            auto reportsppow = std::make_unique<dat::ReportComplexData1D<double,double>>("SQ_SR",m_mirror_kr,m_mirror_usp);
//            reportsppow->setPower(true);
//            PropagatorCV::addReport(std::move(reportsppow));
        }
        SpidaRCVT& m_spi;
        std::vector<dcmplx> m_usp;
        std::vector<dcmplx> m_uphys;
        std::vector<double> m_mirror_r;
        std::vector<double> m_mirror_kr;
        std::vector<double> m_shifted_omega;

        std::vector<dcmplx> m_mirror_uphys;
        std::vector<dcmplx> m_mirror_shift_usp;
};

int main()
{
    unsigned nr = 80;
    double rmax = 4.0;
    BesselRootGridR gridR(nr,rmax);
    unsigned nt = 512;
    double a = -6.0;
    double b = 6.0;
    UniformGridCVT gridT(nt,a,b);

    unsigned num_threads = 4;
    NLSRT model(gridR,gridT,num_threads);

    std::filesystem::path dirpath("nlsRT_propagator_files");
    Propagator propagator(dirpath,model);
    propagator.setStepsPerOutput(4);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(12);

    bool use_refs = true; // Use references to model.L and model.NL, rather than copying
    ETD35 solver(model.L(),model.NL(),use_refs);
    solver.setEpsRel(1e-4);
    solver.setLogProgress(true);
    solver.setLogFrequency(12);
    solver.setNumThreads(num_threads);
    double z0 = 0.0;
    double zf = 0.3;
    solver.evolve(propagator,z0,zf,zf/100.0);
    propagator.report(zf);
    return 0;
}





