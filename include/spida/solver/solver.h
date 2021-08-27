
#ifndef SOLVER_H_
#define SOLVER_H_


#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <memory>
#include <exception>
#include <filesystem>
#include <pwutils/pwstats.h>
#include <pwutils/pwmath.hpp>
#include <pwutils/pwthreads.h>
#include "spida/helper/constants.h"


namespace spida{

const int MAX_LOOP = 100;
using pw::StatCenter;

class ModelCV;
class PropagatorCV;
class ReportCenter;

class SolverCV
{
  public:
      SolverCV(ModelCV& model,PropagatorCV& prs);
      virtual ~SolverCV();
      void setStatFrequency(int val) {m_stat.setReportFrequency(val);}
      void setLogProgress(bool val); 
      void setCurrentTime(double t) {m_tcurrent = t;}
      void setTargetDirectory(const std::filesystem::path& dirpath);

      ModelCV& model() {return m_model;}
      PropagatorCV& propagator() {return m_prs;}
      ReportCenter& reportCenter();
      int size() const;

      StatCenter& statCenter() {return m_stat;}
      int size() {return m_sz;}
      double currentTime() {return m_tcurrent;}
      double dtLast() {return m_dt_last;}
      bool logProgress() {return m_log_progress;}

      virtual void evolve(double t0,double tf,double& dt) = 0;
      void computeCo(double dt);

      void fileReportStats();
      void reportStats();

  private:
      virtual void updateCoefficients([[maybe_unused]] double dt) {};

      ModelCV& m_model;
      PropagatorCV& m_prs;
      std::unique_ptr<ReportCenter> m_report_center;
      StatCenter m_stat;

      double m_tcurrent;
      double m_dt_last;
      bool m_log_progress;
      int m_sz;
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
      double computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew);
      void checkLoopCount(int num_loops);
      void checkStepSize(double step_size);

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
//      double norm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew); // returns s-val, not error
//      double norm1(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew); // returns s-val, not error
//      double weightedNorm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew); // returns s-val, not error
//      double normInf(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew); // returns s-val, not error
//      double normSys(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew); // returns s-val, not error
//      void worker_norm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,long double* esum,long double* ysum);
//      void worker_norm1(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double* esum,double* ysum);
//      void worker_weighted_norm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double,double* esum,double* ysum);
//      void worker_norminf(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double* emax,double* ymax);
//      void worker_normsys(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double cutoff,double* sval);
};


class SolverCV_AS : public SolverCV
{
  public:
      SolverCV_AS(ModelCV& cmodel,PropagatorCV& cprs,double sf,double qv);
      virtual ~SolverCV_AS(); 
  
      void setIncrementThreshold(double val); 
      void setDecrementThreshold(double val);
      void setEpsRel(double val);
      void setNorm(std::string);
      void setAccept(bool val) {m_accept = val;}

      bool accept() {return m_accept;}
  
      void evolve(double t0,double tf,double& h_next);
      void step(double& h,double& h_next);
      
      std::vector<dcmplx>& getY() {return m_yv;}
      std::vector<dcmplx>& getErr() {return m_errv;}

  private:
      virtual void updateCoefficients([[maybe_unused]] double dt) {};
      virtual void updateStages(const std::vector<dcmplx>& in,std::vector<dcmplx>& y,std::vector<dcmplx>& err) = 0;
      std::unique_ptr<Control> m_control;
      std::vector<dcmplx> m_yv;
      std::vector<dcmplx> m_errv;
      bool m_accept;
};


class SolverCV_CS : public SolverCV 
{
  public:
      SolverCV_CS(ModelCV& cmodel,PropagatorCV& cprs);
      virtual ~SolverCV_CS() {};
      void setCountTime(bool val) {m_count_time = val;}

      void step(double h);
      void evolve(double t0,double tf,double& dt);

  private:
      virtual void updateCoefficients([[maybe_unused]] double dt) {};
      virtual void updateStages(std::vector<dcmplx>& in) = 0;
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

#endif


