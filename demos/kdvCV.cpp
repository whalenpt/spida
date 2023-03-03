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

#include <spida/CVX.h>
#include <spida/grid/uniformCVX.h>
#include <spida/helper/constants.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/ETDAS.h>
#include <pwutils/report/dat.hpp>
#include <fstream>

//------------------------------------------------------------------------------

using namespace spida;

/// 
/// @brief Class for holding a linear operator and nonlinear function
/// of the KdV equation for modeling a complex-valued propagating wave.   
/// Also holds SpidaCVX object which has helper functions for ffts and derivatives.
/// 

class KdV
{
    public: 

		/// 
		/// @brief Class constructor copies numerical grid, initializes a SpidaCVX object,
		/// and defines both the linear operator and nonlinear function 
		/// needed by the rkstiff solver. 
		/// @param grid UniformGridCVX object for describing a complex-valued uniform numerical grid 
		/// 

        KdV(const UniformGridCVX& grid) : 
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

        std::vector<dcmplx> L; /**< Linear operator, dcmplx is equivalent to std::complex<double> */

		/// Nonlinear function definition
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
        UniformGridCVX m_grid; /**< Uniformly spaced numerical grid for complex-valued functions */
        SpidaCVX m_spida; /**< FFTs and differentiation functions */
        std::vector<dcmplx> m_uphys; 
        std::vector<dcmplx> m_uxphys;
        std::vector<dcmplx> m_uxsp;
};

//
// Helper class for reporting files based on data generated from the Solver used
//

class PropagatorKdV : public PropagatorCV
{
    public:

		/// 
		/// @brief Constructor creates vectors used during propagation and stores
		//  the propagation model. Determines how the Solver class reports data files
		//  during propagation
		/// @param path Where to output data files generated during propagation
		/// @param md KdV model which contains a SpidaCVX object with fft functions
		
        PropagatorKdV(const std::filesystem::path& path,KdV& md) : 
            PropagatorCV(path),
            m_spi(md.spida()),
            m_usp(md.spida().getGridX().getNsx(),0.0),
            m_uphys(md.spida().getGridX().getNx(),0.0),
			m_shifted_kx(m_spi.getGridX().freqshift(m_spi.getGridX().getSX())),
            m_shifted_usp(m_spi.getGridX().getNsx())
        {

			 // set up initial wave profile (5-solitons)
			 
			 std::vector<double> A0{0.6,0.5,0.4,0.3,0.2};
			 std::vector<double> x0{-120.0,-90.0,-60.0,-30.0,0.0};
			 const std::vector<double>& x  = m_spi.getX();

			 for(auto i = 0; i < x.size(); i++){
				for(auto k = 0; k < A0.size(); k++){
					m_uphys[i] += 0.5*pow(A0[k],2)/pow(cosh(A0[k]*(x[i]-x0[k])/2.0),2);
				}
			 }

			 // convert initial spatial wave profile to spectral space, store in m_usp
			 m_spi.X_To_SX(m_uphys,m_usp);

             // freq shift of spectral components 
			 m_spi.getGridX().freqshift(m_usp,m_shifted_usp); 

			 initReport();
		}

		/// Destructor
        ~PropagatorKdV() {}  
		
        /// @brief Pure virtual function of PropagatorCV that must be implemented.
        /// This function is called before each Solver report (allows for updating of real space fields)
		/// @param t Current propagation time
		
        void updateFields(double t) { 
			m_spi.SX_To_X(m_usp,m_uphys);
			m_spi.getGridX().freqshift(m_usp,m_shifted_usp);
		}

		/// Returns propagating field array
        std::vector<dcmplx>& propagator() {return m_usp;} 

    private:
        /// @brief Helper function that feeds PropagatorCV information on what to report
        void initReport() {

			// Quick access to x vector
            const std::vector<double>& x  = m_spi.getGridX().getX();

            // add report for real-space field
            PropagatorCV::addReport(std::make_unique<dat::ReportComplexData1D<double,double>>("X",x,m_uphys));

            // add report for spectral-space (the propagator)
            PropagatorCV::addReport(std::make_unique<dat::ReportComplexData1D<double,double>>("SX",m_shifted_kx,m_shifted_usp));
        }

        SpidaCVX& m_spi;
        std::vector<dcmplx> m_usp;
        std::vector<dcmplx> m_uphys;
		std::vector<double> m_shifted_kx;
		std::vector<dcmplx> m_shifted_usp;
};


int main()
{
    unsigned nx = 512;
    double minx = -150.0;
    double maxx = 150.0;
    std::string outdir("kdv_files_CV");

    UniformGridCVX grid(nx,minx,maxx);
    KdV model(grid);

    std::filesystem::path dirpath("kdv_files_CV");
    PropagatorKdV propagator(dirpath,model);

    propagator.setStepsPerOutput(16);
    propagator.setLogProgress(true);
    propagator.setLogFrequency(16);

    ETD34 solver(model.L,model.NL);
    solver.setEpsRel(1e-4);
    solver.setLogProgress(true);
    solver.setLogFrequency(16);

    double t0 = 0.0;
    double tf = 600.0;
    double h_init = 0.1;
    solver.evolve(propagator,t0,tf,h_init);

    return 0;

}





