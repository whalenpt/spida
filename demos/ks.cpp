/**------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Complete
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

/// 
/// @brief Class for holding a linear operator and nonlinear function
/// of the Kuramoto–Sivashinsky (KS) equation for modeling a real-valued propagating field.   
/// KS_RV also holds a SpidaRVX object which has helper functions for ffts and derivatives.
/// 

class KS_RV
{
    public: 

		/// 
		/// @brief Class constructor copies numerical grid, initializes a SpidaRVX object,
		/// and defines both the linear operator and nonlinear function 
		/// needed by the rkstiff solver. 
		/// @param grid UniformGridRVX object for describing a real-valued uniform numerical grid 
		/// 

        explicit KS_RV(const UniformGridRVX& grid) : 
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
                    m_spi.SX_To_X(in,m_uphys); // Transform from spectral-space to physical-space
                    m_spi.dSX(in,m_uxsp); // Take derivative of function in spectral-space
                    m_spi.SX_To_X(m_uxsp,m_uxphys); // Convert spectral-space derivative to real-space derivative
                    for(auto i = 0; i < m_grid.getNx(); i++)
                        m_uphys[i] = -m_uphys[i]*m_uxphys[i];
                    m_spi.X_To_SX(m_uphys,out);
                };
            }

        /// Linear operator access 
        std::vector<dcmplx>& L() {return m_L;} 

		/// Nonlinear function access
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)>& NL() {return m_NL;}

        /// Spida object access
        SpidaRVX& spida() {return m_spi;}

    private:
        UniformGridRVX m_grid; /**< Uniformly spaced numerical grid for real-valued functions */
        SpidaRVX m_spi; /**< FFTs and differentiation functions */
        std::vector<double> m_uphys; 
        std::vector<double> m_uxphys;
        std::vector<dcmplx> m_uxsp;
        std::vector<dcmplx> m_L; /**< Linear operator */
        std::function<void(const std::vector<dcmplx>&,std::vector<dcmplx>&)> m_NL; /**< Nonlinear function */
};

///
/// Helper class for reporting files based on data generated from the Solver 
///

class PropagatorKS : public PropagatorCV
{
    public:

		/// 
		/// @brief Constructor creates vectors used during propagation and stores
		//  the propagation model. Determines how the Solver class reports data files
		//  during propagation
		/// @param path Where to output data files generated during propagation
		/// @param md KS_RV model contains a SpidaRVX object with fft functions

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

		/// Destructor
        ~PropagatorKS() override = default;

        /// @brief Pure virtual function of PropagatorCV that must be implemented.
        /// This function is called before each Solver report (allows for updating of real space fields)
		/// @param t Current propagation time

        void updateFields(double t) override { m_spi.SX_To_X(m_usp,m_uphys);}

		/// Returns propagating field array
        std::vector<dcmplx>& propagator() override {return m_usp;}

    private:

        /// @brief Helper function that initializes PropagatorCV with information on what to report
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
    unsigned N = 8192; // Number of grid points
    double a = 0.0; // Numerical grid left end point
    double b = 32.0*spida::PI; // Numerical grid right end point

    UniformGridRVX grid(N,a,b); // Numerical grid
    KS_RV model(grid); // KS model

    std::filesystem::path dirpath("ks_propagator_files"); 
    PropagatorKS propagator(dirpath,model); // Propagator used in solver
    propagator.setStepsPerOutput(5); 
    propagator.setLogProgress(true);
    propagator.setLogFrequency(200);

    ETD34 solver(model.L(),model.NL());
    solver.setEpsRel(1e-4);
    solver.setLogProgress(true);
    solver.setLogFrequency(200);
    solver.evolve(propagator,0.0,50.0,0.5); // Evolve propagator class

    return 0;
}