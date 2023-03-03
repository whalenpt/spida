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
        virtual ~SolverCV() = default;

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
        const LinOp& L() const {return *m_L;}
        /// Accessor for nonlinear function
        const NLfunc& NL() const {return *m_NL;}

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
        double currentTime() const {return m_tcurrent;}

        /// Get the last step size
        double dtLast() const {return m_dt_last;}

        /// Accesses whether logging is enabled
        bool logProgress() const {return m_log_progress;}

        /// File report statistics
        void fileReportStats(const std::filesystem::path& dirpath) const;

        /// Non-file report statistics
        void reportStats() const;

        /// Accessor of the StatCenter
        StatCenter& statCenter() {return m_stat;}

  private:

      ///
      /// @brief Updates coefficients used in the numerical solver
      /// @param dt Current step size of the numerical solver
      ///

      virtual void updateCoefficients([[maybe_unused]] double dt) noexcept {};

      std::unique_ptr<LinOp> m_Lptr; /**< Linear operator */
      std::unique_ptr<NLfunc> m_NLptr; /**< Nonlinear function */
      const LinOp* m_L;
      const NLfunc* m_NL;

      StatCenter m_stat;
      double m_tcurrent{0.0}; /**< Current time */
      double m_dt_last{0.0}; /**< Previous step size */
      bool m_log_progress{false}; /**< Determines whether to log data from propagation */
      bool m_use_refs; /**< Use L and NL passed into class constructor (dont copy) */
      pw::ThreadManager m_thmgt{1}; /**< Holds number of threads and helper functions */
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

        enum class ErrorNorm{NORM1,NORM2,NORMINF,NORMSYS,NORMW2};
        Control(double safetyF,double qv,double epsR,double inF,double decF) :
                    m_safeFact(safetyF),
                    m_q(qv),
                    m_epsRel(epsR),
                    m_incrFact(inF),
                    m_decrFact(decF) {}

        ~Control() = default;
        void setIncrementThreshold(double val); 
        void setDecrementThreshold(double val);
        void setEpsRel(double val);
        void setNorm(ErrorNorm norm) {m_normType = norm;}
        double computeS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) const noexcept;
        double computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) const noexcept;
        bool checkLoopCount(unsigned num_loops) const noexcept;
        bool checkStepSize(double step_size) const noexcept;

    private:
        double m_safeFact;
        double m_q;
        double m_epsRel;
        double m_incrFact;
        double m_decrFact;
        ErrorNorm m_normType{ErrorNorm::NORM2};
};


///  Virtual base class for adaptive-step rkstiff solvers propagating complex-valued fields

class SolverCV_AS : public SolverCV
{

    public:

        /// Initializes an adaptive-step solver class
        SolverCV_AS(const LinOp& L,const NLfunc& NL,double sf,double qv,bool use_refs=false);

        /// Specify descructor in subclasses
        ~SolverCV_AS() override = default; 

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

        bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept override;

        /// @brief Propagates a PropagatorCV object from t0 to tf
        /// @param propagator Class holding the field being propagated plus functions and other arrays
        /// @param t0 Initial time
        /// @param tf Final time
        /// @param h_init Initial time step-size
        /// @return Boolean that says whether the propagation was successful

        bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept override;

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
        /// @param Control::ErrorNorm enum value for norm 
        void setNorm(Control::ErrorNorm norm) {m_control->setNorm(norm);}

        void setAccept(bool val) {m_accept = val;}
        bool accept() const {return m_accept;}

        /// Accessor for the computed updated field
        std::vector<dcmplx>& getY() {return m_yv;}
        /// Accessor for the computed stage error
        std::vector<dcmplx>& getErr() {return m_errv;}

    private:
        void updateCoefficients([[maybe_unused]] double dt) noexcept override {
            // Intentionally unimplemented...
        };
        virtual void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,std::vector<dcmplx>& err) noexcept = 0;
        std::unique_ptr<Control> m_control;
        std::vector<dcmplx> m_yv;
        std::vector<dcmplx> m_errv;
        bool m_accept{false};
        void logStartComputeTime();
        void logEndComputeTime();
        void logStartStepTime();
        void logEndStepTime();
        void logStepRejected(double s); 
        void logStepSizeIncreased(double s); 
};


class SolverCV_CS : public SolverCV 
{
  public:
      using SolverCV::SolverCV;
      ~SolverCV_CS() override = default;
      // takes one numerical step with a step of size h
      void step(std::vector<dcmplx>& u,double h) noexcept;
      // steps the vector u from time t0 to time tf with step size h
      bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept override;
      // evolve with a propagator
      bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept override;

      void setCountTime(bool val) {m_count_time = val;}
  private:
      virtual void updateStages(std::vector<dcmplx>& in) noexcept = 0;
      bool m_count_time{false};
};


class SolverException : public std::exception
{
    public:
        SolverException() = default;
        ~SolverException() override = default;
};

class StepSizeException : public SolverException
{
    public:
        explicit StepSizeException(double val,double minval) 
        {
           m_step_val = std::to_string(val);
           m_minstep_val = std::to_string(minval);
           m_msg = "SOLVER FAILED! The solver step size of " + 
                m_step_val + " is below the specified minimum of " + m_minstep_val + ".";
        }
        ~StepSizeException() override = default;
        const char* what() const noexcept override {
            return m_msg.c_str();
        }
    private:
   	    std::string m_msg;
        std::string m_step_val;
        std::string m_minstep_val;
};


class LoopException : public SolverException
{
    public:
        explicit LoopException(int val) {
					m_maxloops = std::to_string(val);
					m_msg = "SOLVER FAILED! The adaptive step solver has reduced its"\
                " step size " + m_maxloops + " times without achieving an error below epsRel.";
				}
        ~LoopException() override = default;
        const char* what() const noexcept override{
            return m_msg.c_str();
        }
    private:
        std::string m_maxloops;
        std::string m_msg;
};




}



