/**------------------------------------------------------------------------------
 *   
 *    Author: Patrick Whalen   
 *    Email: whalenpt@gmail.com
 *    Status: Complete
 *    Date: 08/27/21
 *    Description: Implementation of kdv PDE with a Propagator class (automated file reporting)
 *
------------------------------------------------------------------------------*/

#pragma once

// HEADERS, INCLUDES, GLOBAL VARS/DECLARATIONS, ETC. 

#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <memory>
#include <exception>
#include <filesystem>
#include <functional>
#include <pwutils/pwstats.h>
#include <pwutils/pwmath.hpp>
#include <pwutils/pwthreads.h>
#include "spida/helper/constants.h"

//------------------------------------------------------------------------------


// define everything in the spida namespace
namespace spida{

/// Maximum number of step size attempts to reach required error threshold
constexpr auto MAX_LOOP = 100;


using pw::StatCenter;
class PropagatorCV;
using NLfunc = std::function<void(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)>; 
using LinOp = std::vector<dcmplx>;


///  Virtual base class for rkstiff solvers for propagating complex-valued fields

class SolverCV
{
    public:

        ///  @brief Initializes SolverCV base class
        ///  @param L Linear operator
        ///  @param NL Nonlinear function

        SolverCV(const LinOp& L,const NLfunc& NL,bool use_refs=false);
        virtual ~SolverCV();

        /// @brief Propagates a vector u from t0 to tf using the class LinOp and NLfunc
        /// @param u Field being propagated
        /// @param t0 Initial time
        /// @param tf Final time
        /// @param h_init Initial time step-size
        virtual bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h_init) noexcept = 0;

        /// @brief Propagates a PropagatorCV object from t0 to tf using the class LinOp and NLfunc
        /// @param propagator Class holding the field being propagated plus functions and other arrays
        /// @param t0 Initial time
        /// @param tf Final time
        /// @param h_init Initial time step-size

        virtual bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept = 0;

        ///
        /// @brief Updates coefficients used in the numerical solver
        /// @param dt Current step size of the numerical solver
        ///

        void computeCo(double dt) noexcept;

        /// Accessor for linear operator
        const LinOp& L() {return *m_L;}
        /// Accessor for nonlinear function
        const NLfunc& NL() {return *m_NL;}

        /// Size of the array of the field being propagated
        unsigned size() const;

        /// Sets logging out to console
        void setLogProgress(bool val) { m_log_progress = val; }

        /// Sets logging frequency
        void setLogFrequency(unsigned val) { m_stat.setLogFrequency(val);};

        /// Sets current time of propagator
        void setCurrentTime(double t) {m_tcurrent = t;}

        /// Sets number of threads that can be utilized by Solver
        void setNumThreads(unsigned val) {m_thmgt.setNumThreads(val);}

        /// Accesses the number of threads that can be utilized by the Solver
        unsigned numThreads() const {return m_thmgt.getNumThreads();}

        /// Accesses the threadManager
        pw::ThreadManager& threadManager() {return m_thmgt;}

        /// Get the current time
        double currentTime() {return m_tcurrent;}

        /// Get the last step size
        double dtLast() {return m_dt_last;}

        /// Accesses whether logging is enabled
        bool logProgress() {return m_log_progress;}

        /// File report statistics
        void fileReportStats(std::filesystem::path& dirpath);

        /// Non-file report statistics
        void reportStats();

        /// Accessor of the StatCenter
        StatCenter& statCenter() {return m_stat;}

  private:

      ///
      /// @brief Updates coefficients used in the numerical solver
      /// @param dt Current step size of the numerical solver
      ///

      virtual void updateCoefficients([[maybe_unused]] double dt) noexcept {};

      const LinOp* m_L; /**< Linear operator */
      const NLfunc* m_NL; /**< Nonlinear function */

      StatCenter m_stat;
      double m_tcurrent; /**< Current time */
      double m_dt_last; /**< Previous step size */
      bool m_log_progress; /**< Determines whether to log data from propagation */
      bool m_use_refs; /**< Use L and NL passed into class constructor (dont copy) */
      pw::ThreadManager m_thmgt; /**< Holds number of threads and helper functions */
};

///  Helper class for computing step updates
class Control{
    public:
        static constexpr double MAX_S = 4.0;
        static constexpr double MIN_S = 0.25;
        static const unsigned MAX_LOOP = 100;
        static constexpr double MIN_H = 1.0e-15;

        /// 
        ///  @brief Constructor for Control class
        ///  @param safetyF Safety factor for predicting best next step size
        ///  @param qv Order coefficient 
        ///  @param epsR Relative error tolerance
        ///  @param inF Increment factor
        ///  @param decF Decrement factor
        ///  @param dim Dimension of system
        ///

        Control(double safetyF,double qv,double epsR,double inF,double decF,int dim);
        ~Control() {}
        void setIncrementThreshold(double val); 
        void setDecrementThreshold(double val);
        void setEpsRel(double val);
        void setNorm(std::string);
        double computeS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew);
        double computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) noexcept;
        bool checkLoopCount(unsigned num_loops) noexcept;
        bool checkStepSize(double step_size) noexcept;

    private:
        double m_safeFact;
        double m_q;
        double m_epsRel;
        double m_incrFact;
        double m_decrFact;
        int m_sz;
        int m_normType;
        enum {NORM1,NORM2,NORMINF,NORMSYS,NORMW2};
};


///  Virtual base class for adaptive-step rkstiff solvers propagating complex-valued fields

class SolverCV_AS : public SolverCV
{

    public:

        /// Initializes an adaptive-step solver class
        SolverCV_AS(const LinOp& L,const NLfunc& NL,double sf,double qv,bool use_refs=false);

        /// Specify descructor in subclasses
        virtual ~SolverCV_AS(); 

        /// @brief Take one numerical step of the vector u with recommended step size h,
        /// but only if h satisfies the error control.
        /// @param u Field to step forward
        /// @param h Step size to take
        /// @param h_next Suggested next step size
        /// @return Boolean determining whether the attempted step was successful

        bool step(std::vector<dcmplx>& u,double& h,double& h_next) noexcept;

        /// @brief Propagates a field u from t0 to tf
        /// @param u Field being propagated 
        /// @param t0 Initial time
        /// @param tf Final time
        /// @param h_init Initial time step-size
        /// @return Boolean that says whether the propagation was successful

        bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept;

        /// @brief Propagates a PropagatorCV object from t0 to tf
        /// @param propagator Class holding the field being propagated plus functions and other arrays
        /// @param t0 Initial time
        /// @param tf Final time
        /// @param h_init Initial time step-size
        /// @return Boolean that says whether the propagation was successful

        bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept;

        /// @brief Setter for increment threshold which dictates how often to increase the solver step size
        /// @param val Increment threshold is at least 1 with a higher value indicating
        ///            less frequent step size increases and vice versa (default is around 1.25)
        
        void setIncrementThreshold(double val); 

        /// @brief Setter for decrement threshold which dictates how often to decrease the solver step size
        /// @param val Decrement threshold is below 1, the step size is reduced by at least this value,
        ///            even if the suggested step size reduction is higher (similar to the safety factor)        

        void setDecrementThreshold(double val);

        /// @brief Setter for relative error tolerance
        /// @param val Relative error tolerance
        
        void setEpsRel(double val);


        /// Setter for error norm, default is the 2-norm
        /// @param str String value of norm 
        /// Note: May want to swap out string with enum class
        void setNorm(std::string str);

        void setAccept(bool val) {m_accept = val;}
        bool accept() {return m_accept;}

        /// Accessor for the computed updated field
        std::vector<dcmplx>& getY() {return m_yv;}
        /// Accessor for the computed stage error
        std::vector<dcmplx>& getErr() {return m_errv;}

    private:
        virtual void updateCoefficients([[maybe_unused]] double dt) noexcept {};
        virtual void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,std::vector<dcmplx>& err) noexcept = 0;
        std::unique_ptr<Control> m_control;
        std::vector<dcmplx> m_yv;
        std::vector<dcmplx> m_errv;
        bool m_accept;
};


class SolverCV_CS : public SolverCV 
{
  public:
      SolverCV_CS(const LinOp& L,const NLfunc& NL,bool use_refs=false);
      virtual ~SolverCV_CS() {};
      // takes one numerical step with a step of size h
      void step(std::vector<dcmplx>& u,double h) noexcept;
      // steps the vector u from time t0 to time tf with step size h
      bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept;
      // evolve with a propagator
      bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept;

      void setCountTime(bool val) {m_count_time = val;}
  private:
      virtual void updateCoefficients([[maybe_unused]] double dt) noexcept {};
      virtual void updateStages(std::vector<dcmplx>& in) noexcept = 0;
      bool m_count_time;
};


class SolverException : public std::exception
{
    public:
        SolverException() {}
        virtual ~SolverException() {};
        virtual const char* what() const noexcept = 0;
};

class StepSizeException : public SolverException
{
    public:
        explicit StepSizeException(double val,double minval) 
        {
           step_val = std::to_string(val);
           minstep_val = std::to_string(minval);
           std::string msg = "SOLVER FAILED! The solver step size of " + 
                step_val + " is below the specified minimum of " + minstep_val + ".";
        }
        ~StepSizeException() {};
        const char* what() const noexcept override {
            return msg.c_str();
        }
    private:
				std::string msg;
        std::string step_val;
        std::string minstep_val;
};


class LoopException : public SolverException
{
    public:
        explicit LoopException(int val) {
					maxloops = std::to_string(val);
					msg = "SOLVER FAILED! The adaptive step solver has reduced its"\
                " step size " + maxloops + " times without achieving an error below epsRel.";
				}
        ~LoopException() {};
        const char* what() const noexcept override{
            return msg.c_str();
        }
    private:
        std::string maxloops;
				std::string msg;
};




}



