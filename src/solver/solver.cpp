
#include "spida/solver/solver.h"
#include "spida/constants.h"
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

SolverDC::SolverDC(ModelDC& model,PropagatorsDC& prs)
 :  m_model(model), m_prs(prs)
{ 
    m_tcurrent = 0.0;
    m_dt_last = 0.0;
    m_log_progress = false;

    m_sz = m_model.specDimSize();

    m_stat.addTimer("Computer Time");
    m_stat.addTimer("Solver Stepping");
    m_stat.addTimer("Timer Coefficient");
    m_stat.addClock("Clock Computer Time");
    m_stat.addClock("Clock Stepping");
    m_stat.addClock("Clock Coefficient");
    m_stat.addCounter("Coefficient Evals",1);
}


ReportCenter& SolverDC::reportCenter()
{
    if(!m_report_center)
        throw pw::Exception("File reporting has not been activated."\
        " To enable file reporting use the SolverDC::setFileReport(std::string directory)"\
        " function. Then specify reporting details through the ReportCenter object"\
        " returned from the function SolverDC::reportCenter().");
    else
        return *m_report_center;
}

void SolverDC::computeCo(double dt)
{
    if(fabs(dt-m_dt_last) > NEAR_ZERO) {
        m_dt_last = dt;
        m_stat.startClock("Clock Coefficient");
        m_stat.startTimer("Timer Coefficient");
        updateCoefficients(dt);
        m_stat.endTimer("Timer Coefficient");
        m_stat.endClock("Clock Coefficient");
        m_stat.incrementCounter("Coefficient Evals");
    }
}

void SolverDC::fileReportStats() {
	std::filesystem::path local_path("solver_stat.dat");
    std::filesystem::path dir_path = reportCenter().dirPath();
    std::filesystem::path full_path = dir_path / local_path;
    std::ofstream fout(full_path.c_str());
    m_stat.report(fout);
    fout.close();
}

void SolverDC::reportStats()
{
    m_stat.report(std::cout);
}

void SolverDC::setTargetDirectory(const std::string& dirname)
{
    if(m_report_center){
        m_report_center.reset();
    }
    if(m_model.numDimensions() == 1){
        m_report_center = std::unique_ptr<ReportCenter1D>(new ReportCenter1D(dirname,1,10000));
    } else if(m_model.numDimensions() == 2) {
        m_report_center = std::unique_ptr<ReportCenter2D>(new ReportCenter2D(dirname,1,1,250));
    } else{
        throw pw::Exception("SolverDC::setFileReport","SolverAS can only handle"\
                "1 or 2 dimensions. ");
    }
    m_report_center->addPropagator(&m_prs);
    m_report_center->setLogProgress(m_log_progress);
}

void SolverDC::setLogProgress(bool val) {
    m_log_progress = val;
    if(m_report_center)
        m_report_center->setLogProgress(val);
}

SolverAS_DC::SolverAS_DC(ModelDC& model,PropagatorsDC& prs,double sf,double qv)
 :  SolverDC(model,prs), m_yv(SolverDC::size()), m_errv(SolverDC::size())
{
    m_accept = false;
    m_control = std::unique_ptr<Control>(new Control(sf,qv,1.0e-3,1.25,0.85,SolverDC::size(),\
                SolverDC::model().threadManager()));
    SolverDC::statCenter().addCounter("Step Size Reductions",1);
    SolverDC::statCenter().addCounter("Step Size Increases",1);
}

SolverAS_DC::~SolverAS_DC()
{
}

void SolverAS_DC::setIncrementThreshold(double val) 
{m_control->setIncrementThreshold(val);}
void SolverAS_DC::setDecrementThreshold(double val)
{m_control->setDecrementThreshold(val);}
void SolverAS_DC::setNorm(std::string str) 
{m_control->setNorm(str);}
void SolverAS_DC::setEpsRel(double val) 
{m_control->setEpsRel(val);}


void SolverAS_DC::step(double& h,double& hNext)
{
    std::vector<dcmplx>& in = SolverDC::propagators().spec();
    statCenter().startClock("Clock Computer Time");
    statCenter().startTimer("Computer Time");
    double h_cur = h;
    m_control->checkStepSize(h_cur);
    double s = 1.0;
    int num_loops = 0;
    bool bigError = true;
    while(bigError){
        num_loops++;
        m_control->checkLoopCount(num_loops);
        computeCo(h_cur);
        statCenter().startTimer("Solver Stepping");
        statCenter().startClock("Clock Stepping");
        updateStages(in,m_yv,m_errv);
        statCenter().endClock("Clock Stepping");
        statCenter().endTimer("Solver Stepping");
        statCenter().incrementCounter("Nonlinear Function Evaluations");
  
        s = m_control->computeS(m_errv,m_yv);
        if(s < 1.0){
            if(SolverDC::logProgress())
                std::cout << "Step rejected with s = " << s << std::endl;
            h_cur = s*h_cur;
            m_accept = false;
            statCenter().incrementCounter("Step Size Reductions");
        }
        else{
            SolverDC::setCurrentTime(SolverDC::currentTime()+h_cur);
            for(int i = 0; i < SolverDC::size(); i++)
                in[i] = m_yv[i];
            if(s > 1.0) {
                hNext = s*h_cur;
                statCenter().incrementCounter("Step Size Increases");
                if(SolverDC::logProgress())
                    std::cout << "step increased with s = " << s << std::endl;
            }
            else
                hNext = h_cur;
            h = h_cur;
            bigError = false;
            m_accept = true;
        }
        m_control->checkStepSize(h_cur);
    }
    statCenter().endClock("Clock Computer Time");
    statCenter().endTimer("Computer Time");
    statCenter().statUpdate();
}

void SolverAS_DC::evolve(double t0,double tf,double& dt)
{
    SolverDC::reportCenter().report(t0);
    SolverDC::setCurrentTime(t0);
    if(dt > (tf-t0))
        dt = tf-t0;

    bool report = true;
    double dt_next;
    while(SolverDC::currentTime() < tf && report){
        step(dt,dt_next);
        report = reportCenter().stepUpdate(SolverDC::currentTime());
        if(SolverDC::currentTime() + dt_next > tf)
            dt = tf - SolverDC::currentTime();
        else
            dt = dt_next;
    }
    reportCenter().report(SolverDC::currentTime());
}

SolverCS_DC::SolverCS_DC(ModelDC& model,PropagatorsDC& prs) :
    SolverDC(model,prs)
{
    m_count_time = true;
}

void SolverCS_DC::step(double dt) 
{
    statCenter().startClock("Clock Computer Time");
    statCenter().startTimer("Computer Time");

    std::vector<dcmplx>& in = SolverDC::propagators().spec();
    computeCo(dt);
    statCenter().startClock("Clock Stepping");
    statCenter().startTimer("Solver Stepping");
    updateStages(in);
    statCenter().endTimer("Solver Stepping");
    statCenter().endClock("Clock Stepping");

    SolverDC::setCurrentTime(SolverDC::currentTime()+dt);
    statCenter().endTimer("Computer Time");
    statCenter().endClock("Clock Computer Time");

    statCenter().incrementCounter("Nonlinear Function Evaluations");
    statCenter().statUpdate();
}

void SolverCS_DC::evolve(double t0,double tf,double& dt)
{
    SolverDC::setCurrentTime(t0);
    reportCenter().report(t0);
    bool report = true;
    while(SolverDC::currentTime() < tf && report){
        if ((SolverDC::currentTime() + dt) > tf)
            dt = tf - SolverDC::currentTime();
        step(dt);
        report = reportCenter().stepUpdate(SolverDC::currentTime());
    }
    reportCenter().report(SolverDC::currentTime());
}


Control::Control(double safetyF,double qv,double epsR,\
        double inF,double decF,int dim,pw::ThreadManager& thmgt)
  : safeFact(safetyF),q(qv),epsRel(epsR),
    incrFact(inF),decrFact(decF),sz(dim), 
    th_manage(thmgt), 
    esum(th_manage.getNumThreads(),0.0), 
    ysum(th_manage.getNumThreads(),0.0)
{
  MAX_S = 2.5;
  MIN_S = 0.4;
  MAX_LOOP = 100;
  MIN_H = 1.0e-15;
  normType = NORM2;
}

void Control::setNorm(std::string str)
{
  if(str == "NORM1" || str == "1NORM" || str == "NORM 1" || str == "1 NORM")
    normType = NORM1;
  else if(str == "NORM2" || str == "2NORM" || str == "2 NORM" || str == "NORM 2")
    normType = NORM2;
  else if(str == "INFNORM" || str == "NORMINF" || str == "INF NORM" || str == "NORM INF")
    normType = NORMINF;
  else if(str == "SYSNORM" || str == "SYS NORM" || str == "NORM SYS" || str == "NORMSYS")
    normType = NORMSYS;
  else if(str == "WEIGHTED NORM 2" || str == "WEIGHTED NORM2" || str == "WEIGHTED 2 NORM" \
      || str == "WEIGHTED NORM2") normType = NORMW2;
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
    incrFact = val;
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
    decrFact = val;
}

void Control::setEpsRel(double val)   
{
    if(val < 0.0){
        std::string msg = "The relative error must be positive."\
             "Please check the epsRel value of " + std::to_string(val); 
        throw pw::Exception("Control::setEpsRel",msg);
    }
    epsRel = val; 
}


double Control::computeS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew)
{
    double s = computeRawS(errVec,ynew);
    if(isinf(s))
        return MIN_S;
    else if(isnan(s))
        return MIN_S;
    double s_stab = safeFact*pow(s,1.0/q);
    if(s_stab < 1.0){
        s_stab = std::max(s_stab,MIN_S);
        if(s_stab > decrFact)
            s_stab = decrFact;
    } else{
        s_stab = std::min(s_stab,MAX_S);
        if(s_stab < incrFact)
            s_stab = 1.0;
    }
    return s_stab;
}

void norm2(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew,\
        int sti,int endi,double* esum_val,double* ysum_val)
{
    *esum_val = 0.0; 
    *ysum_val = 0.0;
    for(int i = sti; i < endi; i++){
        *esum_val += pow(abs(errVec[i]),2);
        *ysum_val += pow(abs(ynew[i]),2);
    }
}

double Control::computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew)
{
    std::fill(esum.begin(),esum.end(),0.0);
    std::fill(ysum.begin(),ysum.end(),0.0);
    unsigned int parts = th_manage.getNumThreads();
    
    std::vector<std::thread*> threads;
    std::vector<unsigned int> bounds = th_manage.getBounds(sz);
    if(normType == NORM2){
        for(unsigned int i = 0; i < parts-1; i++)
            threads.push_back(new std::thread(norm2,std::ref(errVec),std::ref(ynew),bounds[i],bounds[i+1],&esum[i],&ysum[i]));
        norm2(errVec,ynew,bounds[parts-1],bounds[parts],&esum[parts-1],&ysum[parts-1]);
    }
    else{
        throw pw::Exception("Control::computeRawS","Did not recognize error "\
                "control Norm metric (1-norm,2-norm,inf-norm,...)");
    }
    for(auto t : threads)
        t->join();
    for(auto t : threads)
        delete t;
    threads.clear();
        
    double enet = std::accumulate(esum.begin(),esum.end(),0.0);
    double ynet = std::accumulate(ysum.begin(),ysum.end(),0.0);
    if(normType == NORM2){
        enet = sqrt(enet);
        ynet = sqrt(ynet);
    } else{
        throw pw::Exception("Control::computeRawS","Did not recognize error "\
                "control Norm metric (1-norm,2-norm,inf-norm,...)");
    }
    double errTol= epsRel*ynet;
    return errTol/enet;
}

void Control::checkLoopCount(int num_loops)
{
    if(num_loops >= MAX_LOOP){
        throw LoopException(num_loops);
    }
}

void Control::checkStepSize(double step_size)
{
    if(step_size < MIN_H){
        throw StepSizeException(step_size,MIN_H);
    }
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

