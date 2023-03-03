=========
SPIDA
=========

.. _rkstiff: https://github.com/whalenpt/rkstiff
.. _indirect target: rkstiff_

.. _CMake: https://cmake.org
.. _indirect target: CMake_


Spectral integration and differentiation algorithms (SPIDA). This project implements several
high-level wrappers for spectral transforms:

#. Fast-Fourier Transforms (FFTs) 
#. Fast Discrete Cosine Transforms (FDCTs) 
#. Hankel-Transforms  

SPIDA also contains adaptive-step Runge-Kutta numerical propagators for stiff partial-differential-equations (PDEs).
These are the same methods as implemented for the python package rkstiff_. 
Two types of numerical methods are available, integrating factor methods (IF), and exponential time-differencing (ETD).
In particular, the solvers provided are

* ETD35 (5th order ETD with 3rd order embedding)
* ETD34 (4th order ETD with 3rd order embedding) 
* IF34 (4th order IF with 3rd order embedding)
* IF45DP (5th order IF with 4th order embedding)

In general, one should prefer ETD35 as it often has the best speed and stability for diagonal systems.
Non-diagonal methods are not implemented as of yet (see the python package rkstiff_ for this functionality).
A detailed discussion of these solvers is provided in the journal article |article|_.

 .. _article: https://www.sciencedirect.com/science/article/pii/S0021999114006743

 .. |article| replace:: *Exponential time-differencing with embedded Runge–Kutta adaptive step control*

Dependencies
------------

Libraries built with

* `KISS FFT <https://github.com/mborgerding/kissfft>`_
* `Nayuki-Lee DCT <https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms>`_ 
* `pwutils <https://github.com/whalenpt/pwutils>`_
* `Boost <https://www.boost.org>`_ (installed version if available, otherwise a slimmed down boost in the external folder)
* `json11 <https://github.com/dropbox/json11>`_
 
Testing done with

* `GoogleTest <https://github.com/google/googletest>`_

Compile and build

* CMake_

Hankel transform usage
----------------------

Hankel transforms of order zero computed on a Bessel root grid:

.. code-block:: c

    #include <spida/grid/besselR.h>
    #include <spida/shape/shapeR.h>
    #include <spida/transform/hankelR.h>
    #include <pwutils/report/dat.hpp>
    #include <pwutils/report/reporthelper.h>
    #include <iostream>
    #include <fstream>

    int main()
    {
        // dcmplx is short for std::complex<double>
        int nr = 200;

        // Set up a grid of Bessel function roots
        spida::BesselRootGridR gridR(nr,12*w0);
        double A0 = 1.0; double w0 = 1.0;
        // GaussR -> A0*exp(-r^2/w0^2)
        spida::GaussR shapeR(gridR,A0,w0);
        
        // HankelTransformR is a Hankel transform using Bessel function roots as the grid
        spida::HankelTransformR transform(gridR);

        std::vector<double> u = shapeR.shapeRV();
        std::vector<double> v(nr);
        // Transform to spectral space
        transform.R_To_SR(u,v);

        dat::ReportData1D<double,double> in_report("R",gridR.getR(),u);
        dat::ReportData1D<double,double> out_report("SR",gridR.getSR(),v);

        std::ofstream os;
        // can set precision for reporting
        os << std::scientific << std::setprecision(3);
        os << in_report;
        os << std::scientific << std::setprecision(8);
        os << out_report;
        return 0;
    }


The GaussR shape class takes as input: a BesselRootGrid, an amplitude, 
and a width (also, potentially a focus parameter for complex-valued inputs).
The GaussR::shapeRV() function returns real-valued double vector
containing a Gaussian profile.

The HankelTransformR class also takes a BesselRootGrid as an input.
The HankelTransformR::R_To_SR(const std::vector<double>& in, std::vector<double>& out) function
will transform the radial Gaussian profile into spectral space (which also should be Gaussian).

Note that pwutils is used here to output results to files named "R.dat" and "SR.dat".
dat::ReportData1D reports two STL vectors into a two-column data file. The spectral 
results are in the SR.dat file.


rkstiff usage
-------------

.. raw:: html

    <embed>
        Consider the Kuramoto-Sivashinsky (KS) equation: 
        <br>
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
         u<sub>t</sub> = -u<sub>xx</sub> - u<sub>xxxx</sub> - uu<sub>x</sub>. 
         
         Converting to spectral space using a Fourier transform (F) we have 
        <br>
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        v<sub>t</sub> = k<sub>x</sub><sup>2</sup>(1- k<sub>x</sub><sup>2</sup>)v - F \{ F<sup>-1</sup> \{v\} F<sup>-1</sup>\{ i k<sub>x</sub> v\} \} 
        <br>
        where v = F{u}. We can then plug L = k<sub>x</sub><sup>2</sup>(1- k<sub>x</sub><sup>2</sup>), and NL(u) =  - F \{ F<sup>-1</sup> \{v\} F<sup>-1</sup>\{ i k<sub>x</sub> v\} \} into an rkstiff solver and propagate the field u in spectral space, converting back to real space when desired.
        For example, the C++ code may look something like this:
    </embed>
|

.. code-block:: c

    #include <spida/grid/uniformRVX.h>
    #include <spida/RVX.h>
    #include <spida/helper/constants.h>
    #include <spida/rkstiff/ETDAS.h>
    #include <spida/propagator/propagator.h>
    #include <pwutils/report/dat.hpp>
    #include <fstream>

    //------------------------------------------------------------------------------

    using namespace spida;

    // KS model for real-valued physical space fields (spectral space is complex)
    class KS_RV
    {
        public: 
            KS_RV(const UniformGridRVX& grid) : m_grid(grid), m_spi(grid), 
                m_uphys(grid.getNx()), m_uxphys(grid.getNx()), m_uxsp(grid.getNsx()),
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
                PropagatorCV(path), m_spi(md.spida()),
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
            std::vector<dcmplx>& propagator() {return m_usp;}
            // updateFields is a pure virtual function of PropagatorCV and must be implemented 
            // This function is called before each Solver report (allows for updating of real space fields)
            void updateFields(double t) { m_spi.SX_To_X(m_usp,m_uphys);}
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

        std::filesystem::path dirpath("ks_propagator_files");
        PropagatorKS propagator(dirpath,model);
        propagator.setStepsPerOutput(5);
        propagator.setLogProgress(true);
        propagator.setLogFrequency(200);

        ETD34 solver(model.L(),model.NL());
        solver.setEpsRel(1e-4);
        solver.setLogProgress(true);
        solver.setLogFrequency(200);
        solver.evolve(propagator,0.0,50.0,0.5);

        return 0;
    }

The solvers, including ETD34, are instantiated with a diagonal linear operator 
as the first argument (L -> std::vector<std::complex<double>>), 
and a nonlinear function as the second argument (NL -> func(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)).

Here KS_RV is a simple class that holds both the linear and nonlinear operators
along with a SpidaRVX object which contains the real-valued (RV) physical-space
to complex-valued (CV) spectral-space transform on a uniform grid (FFT for real-value fields).
KS_RV also holds several intermediate arrays used in the nonlinear function evaluation.

PropagatorKS is a class that inherits from PropagatorCV which is a container
for a complex-valued propagating field. This class has several helper
functions for convenient file reporting, such the number of steps for the
solver to take before each report and whether to log the solvers progress with
std::cout. In particular, the class has two pure virtual functions

* std::vector<spida::dcmplx>& propagator()
* void updateFields(double t) 

that need to be specified in a subclass. The propagator() function returns
the complex-valued array that is propagated by the solver. The updateFields
function is called right before any file report. Note that none of the solvers
require the use of a PropagatorCV class and can use a std::vector input
directly.

The main function sets up the grid, model, propagator, and solver.
The ETD34 evolve function automatically file reports results based
on the settings provided by the PropagatorCV class.

Installation
------------

If not installed, get and install CMake_.

To build and install the SPIDA library from the terminal use

.. code-block:: none

    git clone https://github.com/whalenpt/spida.git
    cd spida
    cmake -S . -B build
    cd build
    cmake --build . (optional flag --config Release/Debug)
    cmake --install .

Check the `usage <./usage>`_ folder for more information on using the library once installed.

Demos
-----

Check out the demos. These can be built by configuring CMake with
the option SPIDA_DEMOS set to ON. On the command line, in the spida directory,
the configure command is:

.. code-block:: none

    cmake -S . -B build -DSPIDA_DEMOS=ON

Testing
-------

Testing done with GoogleTest. Enable testing by configuring CMake
with the option SPIDA_TEST set to ON. On the command line, in the spida directory,
the configure command is:

.. code-block:: none

    cmake -S . -B build -DSPIDA_TEST=ON

License
-------
This project is licensed under the MIT License - see the `LICENSE <./LICENSE>`_ file for details.

Third-party package dependencies use MIT or similarly permissive licenses

Contact
-------
Patrick Whalen - whalenpt@gmail.com



