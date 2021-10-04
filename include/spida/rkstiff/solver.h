// solver.h
#pragma once

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


namespace spida{

constexpr auto MAX_LOOP = 100;

using pw::StatCenter;

class PropagatorCV;
using NLfunc = std::function<void(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)>; 
using LinOp = std::vector<dcmplx>;

class SolverCV
{
  public:
      SolverCV(const LinOp& L,const NLfunc& NL);
      virtual ~SolverCV();

      // evolve propagates a vector u from t0 to tf using the class LinOp and NLfunc
      // t0  is initial time, tf is final time, h_init is initial step size
      virtual bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h_init) noexcept = 0;

      // evolve with a propagator
      virtual bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept = 0;

      void computeCo(double dt) noexcept;
      LinOp& L() {return m_L;}
      NLfunc& NL() {return m_NL;}
      unsigned size() const;

      void setLogProgress(bool val) { m_log_progress = val; }
      void setLogFrequency(unsigned val) { m_stat.setLogFrequency(val);};
      void setCurrentTime(double t) {m_tcurrent = t;}

      void setNumThreads(unsigned val) {m_thmgt.setNumThreads(val);}
      unsigned numThreads() const {return m_thmgt.getNumThreads();}
      pw::ThreadManager& threadManager() {return m_thmgt;}

      double currentTime() {return m_tcurrent;}
      double dtLast() {return m_dt_last;}

      bool logProgress() {return m_log_progress;}
      void fileReportStats(std::filesystem::path& dirpath);
      void reportStats();
      StatCenter& statCenter() {return m_stat;}

  private:
      virtual void updateCoefficients([[maybe_unused]] double dt) noexcept {};
      LinOp m_L;
      NLfunc m_NL;

      StatCenter m_stat;
      double m_tcurrent;
      double m_dt_last;
      bool m_log_progress;
      pw::ThreadManager m_thmgt;
};

void norm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,double* esum_val,double* ysum_val);

class Control{
  public:
      Control(double safetyF,double qv,double epsR,double inF,double decF,int dim,pw::ThreadManager& cthmgt);
      ~Control() {}
      void setIncrementThreshold(double val); 
      void setDecrementThreshold(double val);
      void setEpsRel(double val);
      void setNorm(std::string);
      void setNumThreads(int numThreads) {th_manage.setNumThreads(numThreads);}
      double computeS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew);
      double computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) noexcept;
      bool checkLoopCount(unsigned num_loops) noexcept;
      bool checkStepSize(double step_size) noexcept;

  private:
      double safeFact;
      double q;
      double epsRel;
      double incrFact;
      double decrFact;
      int sz;
      int normType;
      double MAX_S;
      double MIN_S;
      int MAX_LOOP;
      double MIN_H;
      pw::ThreadManager& th_manage;
      std::vector<double> esum;
      std::vector<double> ysum;
 
      enum {NORM1,NORM2,NORMINF,NORMSYS,NORMW2};
};


class SolverCV_AS : public SolverCV
{
  public:
      SolverCV_AS(const LinOp& L,const NLfunc& NL,double sf,double qv);
      virtual ~SolverCV_AS(); 

      // step takes one numerical step of the vector u with recommended step size h, 
      // though potentially not h (if h is too large to satisfy the error control).
      // Next step is recommended to be of size h_next
      bool step(std::vector<dcmplx>& u,double& h,double& h_next) noexcept;

      // numerically propagates u from t0 to tf with initial step size h
      bool evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept;
  
      // evolve with a propagator
      bool evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept;

      void setIncrementThreshold(double val); 
      void setDecrementThreshold(double val);
      void setEpsRel(double val);
      void setNorm(std::string);
      void setAccept(bool val) {m_accept = val;}
      bool accept() {return m_accept;}
      std::vector<dcmplx>& getY() {return m_yv;}
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
      SolverCV_CS(const LinOp& L,const NLfunc& NL);
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



