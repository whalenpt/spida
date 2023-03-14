/**------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Development
 *    Date: 08/27/21
 *    Description: Implementation of kdv PDE using real-valued physical space values
 *
------------------------------------------------------------------------------*/

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <spida/RVX.h>
#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>
#include <spida/rkstiff/ETDAS.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using spida::dcmplx;
using spida::ii;

///
/// @brief Class for holding a linear operator and nonlinear function
/// of the KdV equation for modeling a real-valued propagating wave.   
/// Also holds SpidaRVX object which has helper functions for ffts and derivatives.
///

class KdV_RV
{
    public: 

		/// 
		/// Class constructor copies numerical grid, initializes a SpidaRVX object,
		/// and defines both the linear operator and nonlinear function 
		/// needed by the rkstiff solver. 
		/// @param grid UniformGridRVX object for describing a real-valued uniform numerical grid 
		/// 

        explicit KdV_RV(const spida::UniformGridRVX& grid) : 
            L(grid.getNsx()),
            m_grid(grid), 
            m_spida(grid), 
            m_uphys(m_grid.getNx()),
            m_uxphys(m_grid.getNx()),
            m_uxsp(grid.getNsx())
            {
                const auto& sx = grid.getSX();
                for(size_t i = 0; i < sx.size(); i++)
                    L[i] = ii*pow(sx[i],3);
                NL = [this](const std::vector<dcmplx>& in,std::vector<dcmplx>& out){
                    m_spida.SX_To_X(in,m_uphys);
                    m_spida.dSX(in,m_uxsp);
                    m_spida.SX_To_X(m_uxsp,m_uxphys);
                    for(unsigned i = 0; i < m_grid.getNx(); i++)
                        m_uphys[i] = -6.0*m_uphys[i]*m_uxphys[i];
                    m_spida.X_To_SX(m_uphys,out);
                };
            }

        std::vector<dcmplx> L; /**< Linear operator, dcmplx is equivalent to std::complex<double> */
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> NL; /**< Nonlinear function */
        spida::SpidaRVX& spida() {return m_spida;}

        spida::UniformGridRVX m_grid; /**< Grid class holding both real-space uniform grid and spectral-space grid */
        spida::SpidaRVX m_spida; /** < Spida object which contains ffts and differentiation functions */
        std::vector<double> m_uphys; /**< Physical space field */
        std::vector<double> m_uxphys; /**< Spatial derivative */
        std::vector<dcmplx> m_uxsp; /**< Spatial derivative in spectral-space */
};


int main()
{
    unsigned nx = 512; // Number of points in numerical grid
    double minx = -150.0; // Minimum of X-grid. 
    double maxx = 150.0; // Maximum of X-grid. 
    std::filesystem::path outdir("kdv_files_RV"); // Output directory for simulation 

    spida::UniformGridRVX grid(nx,minx,maxx); // Construct numerical grid. 
    KdV_RV model(grid); // Construct KdV_RV model object 
    spida::ETD34 solver(model.L,model.NL); // Setup up ETD34 solver
    solver.setEpsRel(1e-4); // Specify error tolerance 

	// set up initial wave profile (5-solitons)
    std::vector<double> u0(nx,0.0);
    std::vector<double> A0{0.6,0.5,0.4,0.3,0.2};
    std::vector<double> x0{-120.0,-90.0,-60.0,-30.0,0.0};
    const auto& x  = grid.getX();
    for(unsigned i = 0; i < nx; i++){
        for(size_t k = 0; k < A0.size(); k++){
            u0[i] += 0.5*pow(A0[k],2)/pow(cosh(A0[k]*(x[i]-x0[k])/2.0),2);
        }
    }

	// convert initial spatial wave profile to spectral space, store in usp
    std::vector<dcmplx> usp(grid.getNsx());
    model.spida().X_To_SX(u0,usp);

    double t0 = 0.0; // initial time
    double tf = 600.0; // final time
    double h = 0.1; // initial step size
    double h_next = 0.1; // suggested next step size
    double t = t0; // current propagation time

	// copy initial field to uphys
    std::vector<double> uphys(nx);
    std::copy(std::cbegin(u0),std::cend(u0),std::begin(uphys));

	// report initial physical field
    dat::ReportData1D report{"X_0",x,uphys};
    report.setDirPath(outdir);
    report.setItem("t",t0);
    std::cout << "First physical space report file location: " << report.path() << std::endl;

	// report initial spectral field
    dat::ReportData1D reportS{"SX_0",grid.getSX(),usp};
    reportS.setDirPath(outdir);
    std::cout << "First spectral space report file location: " << reportS.path() << std::endl;

    std::ofstream os;
    os << std::scientific << std::setprecision(8);
    os << report;
    os << reportS;

	// initialize counters for number of steps and number of reports
    unsigned step_count = 0;
    unsigned report_count = 1;

	// propagate field with ETD34 class until the propagated time is equal to the final time
    while(t < tf){

        if(!solver.step(usp,h,h_next)){ // propagate unless solver fails
            std::cout << "Step failed." << std::endl;
            return 1;
        }

        t += h; // update current time
        h = h_next; // set next step h to the suggest step h_next

        // Output results every 16th step
        if(step_count % 16 == 0){
			model.spida().SX_To_X(usp,uphys); // Convert spectral-field to physical-field
            report.setName("X_" + std::to_string(report_count));
            reportS.setName("SX_" + std::to_string(report_count));
            report.setItem("t",t);
            reportS.setItem("t",t);
            os << report;
            os << reportS;
            report_count++;
        }
        step_count++;
    }
	model.spida().SX_To_X(usp,uphys); // Convert spectral-field to physical-field
    report.setName("X_" + std::to_string(report_count));
    reportS.setName("SX_" + std::to_string(report_count));
    report.setItem("t",t);
    reportS.setItem("t",t);
    os << report;
    os << reportS;
    os.close();
    return 0;
}