
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

class PropagatorCV;
class ReportCenter;
class ModelCV;

class SolverCV
{
  public:
      SolverCV(ModelCV* model);
      virtual ~SolverCV();

      virtual void evolve(std::vector<dcmplx>& u,double t0,double tf,double& dt) = 0;
      void computeCo(double dt);
      ModelCV& model() {return *m_model;}
      int size() const;

      void setFileReport(std::unique_ptr<PropagatorCV> pr,const std::filesystem::path& dirpath);
      void setTargetDirectory(const std::filesystem::path& dirpath);
      bool fileReportOn() const {return (m_report_center ? true : false);}

      ReportCenter* reportCenter() {return m_report_center.get();}
      void setStatFrequency(int val) {m_stat.setReportFrequency(val);}
      void setLogProgress(bool val); 
      void setCurrentTime(double t) {m_tcurrent = t;}
      double currentTime() {return m_tcurrent;}
      double dtLast() {return m_dt_last;}
      bool logProgress() {return m_log_progress;}

      void fileReportStats();
      void reportStats();
      StatCenter& statCenter() {return m_stat;}

  private:
      virtual void updateCoefficients([[maybe_unused]] double dt) {};
      ModelCV* m_model;
      std::unique_ptr<ReportCenter> m_report_center;
      StatCenter m_stat;
      double m_tcurrent;
      double m_dt_last;
      bool m_log_progress;
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
};


class SolverCV_AS : public SolverCV
{
  public:
      SolverCV_AS(ModelCV* cmodel,double sf,double qv);
      virtual ~SolverCV_AS(); 
      void evolve(std::vector<dcmplx>& u,double t0,double tf,double& h_next);
      void step(std::vector<dcmplx>& u,double& h,double& h_next);
  
      void setIncrementThreshold(double val); 
      void setDecrementThreshold(double val);
      void setEpsRel(double val);
      void setNorm(std::string);
      void setAccept(bool val) {m_accept = val;}
      bool accept() {return m_accept;}
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
      SolverCV_CS(ModelCV* cmodel);
      virtual ~SolverCV_CS() {};
      void step(std::vector<dcmplx>& u,double h);
      void evolve(std::vector<dcmplx>& u,double t0,double tf,double& dt);
      void setCountTime(bool val) {m_count_time = val;}
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


