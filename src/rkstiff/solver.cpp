
#include "spida/rkstiff/solver.h"
#include "spida/helper/constants.h"
#include "spida/propagator/propagator.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <pwutils/pwexcept.h>
#include <pwutils/pwthreads.h>

#include <thread>
#include <numeric>
#include <memory>

namespace spida{

SolverCV::SolverCV(const LinOp& L,const NLfunc& NL,bool use_refs) :
    m_use_refs(use_refs),
    m_thmgt(1)
{ 
    if(use_refs){
        m_L = &L;
        m_NL = &NL;
    } else {
        m_L = new LinOp(L);
        m_NL = new NLfunc(NL);
    }
    m_tcurrent = 0.0;
    m_dt_last = 0.0;
    m_log_progress = false;
    m_stat.setLogFrequency(1);
}

SolverCV::~SolverCV() {
    if(!m_use_refs){
        delete m_L;
        delete m_NL;
    }
}


unsigned SolverCV::size() const
{
    return m_L->size();
}

void SolverCV::computeCo(double dt) noexcept
{
    if(fabs(dt-m_dt_last) > NEAR_ZERO) {
        m_dt_last = dt;
        if(m_log_progress){
            m_stat.startClock("Clock Coefficient");
            m_stat.startTimer("Timer Coefficient");
        }
        updateCoefficients(dt);
        if(m_log_progress){
            m_stat.endTimer("Timer Coefficient");
            m_stat.endClock("Clock Coefficient");
            m_stat.incrementCounter("Coefficient Evals");
        }
    }
}

void SolverCV::fileReportStats(std::filesystem::path& dirpath) {
    std::filesystem::path full_path = dirpath / std::filesystem::path("solver_stat.dat");
    std::ofstream fout(full_path);
    m_stat.report(fout);
    fout.close();
}

void SolverCV::reportStats()
{
    m_stat.report(std::cout);
}


SolverCV_AS::SolverCV_AS(const LinOp& L,const NLfunc& NL,double sf,double qv,bool use_refs)
 :  SolverCV(L,NL,use_refs), m_yv(SolverCV::size()), m_errv(SolverCV::size())
{
    m_accept = false;
    m_control = std::unique_ptr<Control>(new Control(sf,qv,1.0e-3,1.25,0.85,SolverCV::size()));
}

SolverCV_AS::~SolverCV_AS() { }

void SolverCV_AS::setIncrementThreshold(double val) 
{m_control->setIncrementThreshold(val);}
void SolverCV_AS::setDecrementThreshold(double val)
{m_control->setDecrementThreshold(val);}
void SolverCV_AS::setNorm(std::string str) 
{m_control->setNorm(str);}
void SolverCV_AS::setEpsRel(double val) 
{m_control->setEpsRel(val);}


bool SolverCV_AS::step(std::vector<dcmplx>& in,double& h,double& hNext) noexcept
    
{
    if(SolverCV::logProgress()){
        statCenter().startClock("Clock Computer Time");
        statCenter().startTimer("Computer Time");
    }
    double h_cur = h;
    if(!(m_control->checkStepSize(h_cur)))
        return false;
            
    double s = 1.0;
    int num_loops = 0;
    bool bigError = true;
    while(bigError){
        num_loops++;
        if(!(m_control->checkLoopCount(num_loops)))
            return false;
        computeCo(h_cur);
        if(SolverCV::logProgress()){
            statCenter().startTimer("Solver Stepping");
            statCenter().startClock("Clock Stepping");
        }
        updateStages(in,m_yv,m_errv);
        if(SolverCV::logProgress()){
            statCenter().endClock("Clock Stepping");
            statCenter().endTimer("Solver Stepping");
        }
  
        s = m_control->computeS(m_errv,m_yv);
        if(s < 1.0){
            h_cur = s*h_cur;
            m_accept = false;
            if(SolverCV::logProgress()){
                std::cout << "Step rejected with s = " << s << std::endl;
                statCenter().incrementCounter("Step Size Reductions");
            }
        }
        else{
            SolverCV::setCurrentTime(SolverCV::currentTime()+h_cur);
            for(int i = 0; i < SolverCV::size(); i++)
                in[i] = m_yv[i];
            if(s > 1.0) {
                hNext = s*h_cur;
                if(SolverCV::logProgress()){
                    std::cout << "step increased with s = " << s << std::endl;
                    statCenter().incrementCounter("Step Size Increases");
                }
            }
            else
                hNext = h_cur;
            h = h_cur;
            bigError = false;
            m_accept = true;
        }
        if(!(m_control->checkStepSize(h_cur)))
            return false;
    }
    if(SolverCV::logProgress()){
        statCenter().endClock("Clock Computer Time");
        statCenter().endTimer("Computer Time");
        statCenter().statUpdate();
    }
    return true;
}

bool SolverCV_AS::evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept
{
    SolverCV::setCurrentTime(t0);
    if(h > (tf-t0))
        h = tf-t0;
    double dt = h;
    double dt_next;
    bool report = true;
    while(SolverCV::currentTime() < tf && report){
        if(!SolverCV_AS::step(u,dt,dt_next))
            return false;
        if(SolverCV::currentTime() + dt_next > tf)
            dt = tf - SolverCV::currentTime();
        else
            dt = dt_next;
    }
    if(SolverCV::logProgress())
        statCenter().report();
    return true;
}

bool SolverCV_AS::evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept
{
    propagator.report(t0);
    SolverCV::setCurrentTime(t0);
    if(h > (tf-t0))
        h = tf-t0;
    double dt = h;
    double dt_next;
    bool report = true;
    while(SolverCV::currentTime() < tf && report){
        if(!SolverCV_AS::step(propagator.propagator(),dt,dt_next))
            return false;
        report = propagator.stepUpdate(SolverCV::currentTime());
        if(SolverCV::currentTime() + dt_next > tf)
            dt = tf - SolverCV::currentTime();
        else
            dt = dt_next;
    }
    propagator.report(SolverCV::currentTime());
    if(SolverCV::logProgress())
        statCenter().report();
    if(propagator.logProgress())
        propagator.reportStats();
    return true;
}


SolverCV_CS::SolverCV_CS(const LinOp& L,const NLfunc& NL,bool use_refs) :
    SolverCV(L,NL,use_refs)
{
    m_count_time = true;
}

void SolverCV_CS::step(std::vector<dcmplx>& in,double dt) noexcept
{
    if(SolverCV::logProgress()){
        statCenter().startClock("Clock Computer Time");
        statCenter().startTimer("Computer Time");
    }

    computeCo(dt);

    if(SolverCV::logProgress()){
        statCenter().startClock("Clock Stepping");
        statCenter().startTimer("Solver Stepping");
    }

    updateStages(in);
    SolverCV::setCurrentTime(SolverCV::currentTime()+dt);

    if(SolverCV::logProgress()){
        statCenter().endTimer("Solver Stepping");
        statCenter().endClock("Clock Stepping");
        statCenter().endTimer("Computer Time");
        statCenter().endClock("Clock Computer Time");
        statCenter().incrementCounter("Nonlinear Function Evaluations");
        statCenter().statUpdate();

    }
}

bool SolverCV_CS::evolve(PropagatorCV& propagator,double t0,double tf,double h) noexcept
{
    SolverCV::setCurrentTime(t0);
    propagator.report(t0);
    bool report = true;
    while(SolverCV::currentTime() < tf && report){
        if ((SolverCV::currentTime() + h) > tf)
            h = tf - SolverCV::currentTime();
        SolverCV_CS::step(propagator.propagator(),h);
        report = propagator.stepUpdate(SolverCV::currentTime());
    }
    propagator.report(SolverCV::currentTime());
    if(SolverCV::logProgress())
        statCenter().report();
    if(propagator.logProgress())
        propagator.reportStats();
    return true;
}

bool SolverCV_CS::evolve(std::vector<dcmplx>& u,double t0,double tf,double h) noexcept
{
    SolverCV::setCurrentTime(t0);
    bool report = true;
    while(SolverCV::currentTime() < tf && report){
        if ((SolverCV::currentTime() + h) > tf)
            h = tf - SolverCV::currentTime();
        SolverCV_CS::step(u,h);
    }
    if(SolverCV::logProgress())
        statCenter().report();
    return true;
}

Control::Control(double safetyF,double qv,double epsR,\
        double inF,double decF,int dim)
  : m_safeFact(safetyF),m_q(qv),m_epsRel(epsR),
    m_incrFact(inF),m_decrFact(decF),m_sz(dim)
{
  m_normType = NORM2;
}

void Control::setNorm(std::string str)
{
  if(str == "NORM1" || str == "1NORM" || str == "NORM 1" || str == "1 NORM")
    m_normType = NORM1;
  else if(str == "NORM2" || str == "2NORM" || str == "2 NORM" || str == "NORM 2")
    m_normType = NORM2;
  else if(str == "INFNORM" || str == "NORMINF" || str == "INF NORM" || str == "NORM INF")
    m_normType = NORMINF;
  else if(str == "SYSNORM" || str == "SYS NORM" || str == "NORM SYS" || str == "NORMSYS")
    m_normType = NORMSYS;
  else if(str == "WEIGHTED NORM 2" || str == "WEIGHTED NORM2" || str == "WEIGHTED 2 NORM" \
      || str == "WEIGHTED NORM2") m_normType = NORMW2;
  else{
        throw pw::Exception("Control::setNorm","Did not recognize error control Norm string"\
                " (1 norm,2 norm,inf norm,...) ");
  }
}

void Control::setIncrementThreshold(double val)   
{
    if(val > MAX_S){
        std::string msg = "IncrementThreshold must be less than MAX_S = "\
             + std::to_string(MAX_S) + ". Please set the incrementThreshold of "\
             + std::to_string(val) + "to a lower value";
        throw pw::Exception("Control::setIncrementThreshold",msg);
    }
    if(val < 1.0){
        std::string msg = "Failed to set an incrementThreshold of " + std::to_string(val) +\
                " because this parameter must be greater than or equal to 1";
        throw pw::Exception("Control::setIncrementThreshold",msg);
    }
    m_incrFact = val;
}

void Control::setDecrementThreshold(double val)   
{
    if(val < MIN_S){
        std::string msg = "DecrementThreshold must be greater than MIN_S = "\
             + std::to_string(MAX_S) + ". Please set the decrementThreshold of "\
             + std::to_string(val) + "to a higher value";
        throw pw::Exception("Control::setDecrementThreshold",msg);
    }
    if(val >= 1.0){
        std::string msg = "Failed to set a decrementThreshold of " + std::to_string(val) +\
                " because this parameter must be less than 1";
        throw pw::Exception("Control::setDecrementThreshold",msg);
    }
    m_decrFact = val;
}

void Control::setEpsRel(double val)   
{
    if(val < 0.0){
        std::string msg = "The relative error must be positive."\
             "Please check the epsRel value of " + std::to_string(val); 
        throw pw::Exception("Control::setEpsRel",msg);
    }
    m_epsRel = val; 
}


double Control::computeS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew)
{
    double s = computeRawS(errVec,ynew);
    if(std::isinf(s))
        return MIN_S;
    else if(std::isnan(s))
        return MIN_S;
    double s_stab = m_safeFact*pow(s,1.0/m_q);
    if(s_stab < 1.0){
        s_stab = std::max(s_stab,MIN_S);
        if(s_stab > m_decrFact)
            s_stab = m_decrFact;
    } else{
        s_stab = std::min(s_stab,MAX_S);
        if(s_stab < m_incrFact)
            s_stab = 1.0;
    }
    return s_stab;
}

double Control::computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) noexcept
{
    double enorm = 0.0, ynorm = 0.0;
    if(m_normType == NORM2){
        for(size_t i = 0; i < errVec.size(); i++){
            enorm += pow(abs(errVec[i]),2);
            ynorm += pow(abs(ynew[i]),2);
        }
        enorm = sqrt(enorm);
        ynorm = sqrt(ynorm);
    }
    else{
        for(size_t i = 0; i < errVec.size(); i++){
            enorm += pow(abs(errVec[i]),2);
            ynorm += pow(abs(ynew[i]),2);
        }
        enorm = sqrt(enorm);
        ynorm = sqrt(ynorm);
    }

    double errTol= m_epsRel*ynorm;
    return errTol/enorm;
}

bool Control::checkLoopCount(unsigned num_loops) noexcept
{
    if(num_loops >= MAX_LOOP){
        return false;
    }
    return true;
}

bool Control::checkStepSize(double step_size) noexcept
{
    if(step_size < MIN_H){
        return false;
    }
    return true;
}



/*
void Control::worker_norm1(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double* esum,double* ysum) 
{
  double loc_esum = 0.0; 
  double loc_ysum = 0.0;

  for(int i = sti; i < endi+1; i++){
      loc_esum += abs(errVec[i]);
      loc_ysum += abs(ynew[i]);
  }
  *esum = loc_esum;
  *ysum = loc_ysum;
}

void Control::worker_norminf(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double* emax,double* ymax) 
{
  double loc_emax = 0.0; 
  double loc_ymax = 0.0;

  for(int i = sti; i < endi+1; i++){
    double errVal = abs(errVec[i]);
    if(errVal > loc_emax)
      loc_emax = errVal;
    double yval = abs(ynew[i]);
    if(yval > loc_ymax)
      loc_ymax = yval;
  }
  *emax = loc_emax;
  *ymax = loc_ymax;
}


void Control::worker_normsys(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double cutoff,double* sval) 
{
  double loc_epsRel = epsRel; 
  double loc_sval = 1.0e16;
  for(int i = sti; i < endi+1; i++){
    double yval = abs(ynew[i]);
    if(yval > cutoff){
      double temp_s = loc_epsRel*yval/abs(errVec[i]);
      if(temp_s < loc_sval)
        loc_sval = temp_s; 
    }
  }
  *sval = loc_sval;
}

void Control::worker_weighted_norm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,int sti,int endi,int tid,double cutoff,double* esq_sum,double* ysq_sum) 
{
  double loc_esq_sum = 0.0; 
  double loc_ysq_sum = 0.0;
  for(int i = sti; i < endi+1; i++){
    double yval = abs(ynew[i]);
    if(yval > cutoff){
      loc_esq_sum += pow(abs(errVec[i]),2);
      loc_ysq_sum += pow(abs(ynew[i]),2);
    }
  }
  *esq_sum = loc_esq_sum;
  *ysq_sum = loc_ysq_sum;
}

*/



}

