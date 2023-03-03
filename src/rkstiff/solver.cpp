
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
    m_use_refs(use_refs)
{ 
    if(use_refs){
        m_L = &L;
        m_NL = &NL;
    } else {
        m_Lptr = std::make_unique<LinOp>(L);
        m_NLptr = std::make_unique<NLfunc>(NL);
        m_L = m_Lptr.get();
        m_NL = m_NLptr.get();
    }
    m_stat.setLogFrequency(1);
}

unsigned SolverCV::size() const
{
    return static_cast<unsigned>(m_L->size());
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

void SolverCV::fileReportStats(const std::filesystem::path& dirpath) const {
    std::filesystem::path full_path = dirpath / std::filesystem::path("solver_stat.dat");
    std::ofstream fout(full_path);
    m_stat.report(fout);
    fout.close();
}

void SolverCV::reportStats() const
{
    m_stat.report(std::cout);
}

SolverCV_AS::SolverCV_AS(const LinOp& L,const NLfunc& NL,double sf,double qv,bool use_refs)
 :  SolverCV(L,NL,use_refs), 
    m_control(std::make_unique<Control>(sf,qv,1.0e-3,1.25,0.85)),
    m_yv(SolverCV::size()), 
    m_errv(SolverCV::size()) { }

void SolverCV_AS::setIncrementThreshold(double val) {
    m_control->setIncrementThreshold(val);
}

void SolverCV_AS::setDecrementThreshold(double val) {
    m_control->setDecrementThreshold(val);
}

void SolverCV_AS::setEpsRel(double val) {
    m_control->setEpsRel(val);
}

void SolverCV_AS::logStartComputeTime() {
    if(SolverCV::logProgress()){
        SolverCV::statCenter().startClock("Clock Computer Time");
        SolverCV::statCenter().startTimer("Computer Time");
    }
}

void SolverCV_AS::logEndComputeTime() {
    if(SolverCV::logProgress()){
        SolverCV::statCenter().endClock("Clock Computer Time");
        SolverCV::statCenter().endTimer("Computer Time");
    }
}

void SolverCV_AS::logStartStepTime() {
    if(SolverCV::logProgress()){
        SolverCV::statCenter().startTimer("Solver Stepping");
        SolverCV::statCenter().startClock("Clock Stepping");
    }
}

void SolverCV_AS::logEndStepTime() {
    if(SolverCV::logProgress()){
        SolverCV::statCenter().endClock("Clock Stepping");
        SolverCV::statCenter().endTimer("Solver Stepping");
    }
}

void SolverCV_AS::logStepRejected(double s) {
    if(SolverCV::logProgress()){
        std::cout << "Step rejected with s = " << s << std::endl;
        SolverCV::statCenter().incrementCounter("Step Size Reductions");
    }
}

void SolverCV_AS::logStepSizeIncreased(double s) {
    if(SolverCV::logProgress()){
        std::cout << "step increased with s = " << s << std::endl;
        SolverCV::statCenter().incrementCounter("Step Size Increases");
    }
}

bool SolverCV_AS::step(std::vector<dcmplx>& in,double& h,double& hNext) noexcept
{
    this->logStartComputeTime();
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
        this->logStartStepTime();
        updateStages(in,m_yv,m_errv);
        this->logEndStepTime();
        s = m_control->computeS(m_errv,m_yv);
        if(s < 1.0){
            h_cur = s*h_cur;
            m_accept = false;
            this->logStepRejected(s);
        }
        else{
            SolverCV::setCurrentTime(SolverCV::currentTime()+h_cur);
            for(int i = 0; i < SolverCV::size(); i++)
                in[i] = m_yv[i];
            if(s > 1.0) {
                hNext = s*h_cur;
                this->logStepSizeIncreased(s);
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
    this->logEndComputeTime();
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


double Control::computeS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) const noexcept
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

double Control::computeRawS(std::vector<dcmplx>& errVec,std::vector<dcmplx>& ynew) const noexcept
{
    double enorm = 0.0;
    double ynorm = 0.0;
    if(m_normType == ErrorNorm::NORM2){
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

bool Control::checkLoopCount(unsigned num_loops) const noexcept
{
    if(num_loops >= MAX_LOOP){
        return false;
    }
    return true;
}

bool Control::checkStepSize(double step_size) const noexcept 
{
    if(step_size < MIN_H){
        return false;
    }
    return true;
}

}