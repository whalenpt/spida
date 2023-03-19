
#include <algorithm>
#include <vector>
#include <thread>
#include "spida/transform/hankelfftRRVT.h"
#include "spida/transform/fftRVT.h"
#include "spida/transform/hankelR.h"
#include "spida/grid/besselR.h" 
#include "spida/grid/uniformRVT.h" 
#include "spida/helper/funcs.hpp"

namespace spida {

  // Transform using class specific thread pool (not much speed bump, but some)
HankelFFTRRVT::HankelFFTRRVT(const BesselRootGridR& gridR,const UniformGridRVT& gridT,unsigned threads) :
    m_nr(gridR.getNr()),
    m_nt(gridT.getNt()),
    m_nst(gridT.getNst()),
    m_nThreads(threads),
    m_rr(m_nr*m_nt),
    m_sr(m_nr*m_nt),
    m_rs(m_nr*m_nst),
    m_ss(m_nr*m_nst)
{
    for (unsigned m = 1; m < m_nThreads; m++)
    {
        m_ready.push_back(false);
        m_thread.push_back(std::thread(&HankelFFTRRVT::worker_thread,this,m));
    }

    for(unsigned i = 0; i < m_nThreads; i++) {
        m_transformT.push_back(std::make_unique<FFTRVT>(gridT));
        m_transformR.push_back(std::make_unique<HankelTransformR>(gridR));
    }
}

HankelFFTRRVT::~HankelFFTRRVT(){
    setState(State::Done);
    for(auto& thread : m_thread)
        thread.join();
    m_thread.clear();
}

// Updated (m_rr->m_ss)
void HankelFFTRRVT::worker_T_To_ST(unsigned tid)
{
    for (unsigned i = tid; i < m_nr; i+=m_nThreads) 
        m_transformT[tid]->T_To_ST(m_rr.data()+m_nt*i,m_ss.data()+m_nst*i);
}

// Updated (m_ss->m_rr)
void HankelFFTRRVT::worker_ST_To_T(unsigned tid)
{
    for (unsigned i = tid; i < m_nr; i+= m_nThreads)
        m_transformT[tid]->ST_To_T(m_ss.data()+m_nst*i,m_rr.data()+m_nt*i);
}

// New (m_rs->m_ss)
void HankelFFTRRVT::worker_CVT_To_ST(unsigned tid)
{
    for(unsigned i = tid; i < m_nr; i+= m_nThreads)
        m_transformT[tid]->CVT_To_ST(m_rs.data()+m_nt*i,m_ss.data()+m_nst*i);
}

// New (m_ss->m_rs)
void HankelFFTRRVT::worker_ST_To_CVT(unsigned tid)
{
    for (unsigned i = tid; i < m_nr; i+= m_nThreads)
        m_transformT[tid]->ST_To_CVT(m_ss.data()+m_nst*i,m_rs.data()+m_nt*i);
}

// Updated (m_rr->m_sr)
void HankelFFTRRVT::worker_R_To_SR(unsigned tid)
{
    for(unsigned i = tid; i < m_nt; i+=m_nThreads)
        m_transformR[tid]->R_To_SR(m_rr.data()+m_nr*i,m_sr.data()+m_nr*i);
}

// updated (m_sr->m_rr)
void HankelFFTRRVT::worker_SR_To_R(unsigned tid)
{
    for (unsigned j = tid; j < m_nt; j+=m_nThreads)  
        m_transformR[tid]->SR_To_R(m_sr.data()+m_nr*j,m_rr.data()+m_nr*j);
}

// New (m_rs->m_ss)
void HankelFFTRRVT::workerCMP_R_To_SR(unsigned tid)
{
    for(unsigned i = tid; i < m_nst; i+=m_nThreads)
        m_transformR[tid]->R_To_SR(m_rs.data()+m_nr*i,m_ss.data()+m_nr*i);
}

// New (m_ss->m_rs)
void HankelFFTRRVT::workerCMP_SR_To_R(unsigned tid)
{
    for(unsigned i = tid; i < m_nst; i+=m_nThreads)
        m_transformR[tid]->SR_To_R(m_ss.data()+m_nr*i,m_rs.data()+m_nr*i);
}

// Updated
void HankelFFTRRVT::RT_To_SRT(const std::vector<double>& in,std::vector<double>& out) 
{
    transpose(in.data(),m_rr.data(),m_nr,m_nt);
    setState(State::R_To_SR); // (m_rr->m_sr)
    worker_R_To_SR(0);
    setState(State::Wait);
    transpose(m_sr.data(),out.data(),m_nt,m_nr);
}                                  

void HankelFFTRRVT::SRT_To_RT(const std::vector<double>& in,std::vector<double>& out) 
{
    transpose(in.data(),m_sr.data(),m_nr,m_nt);
    setState(State::SR_To_R); // (m_sr->m_rr)
    worker_SR_To_R(0);
    setState(State::Wait);
    transpose(m_rr.data(),out.data(),m_nt,m_nr);
}


void HankelFFTRRVT::RST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out)
{
    std::copy(std::begin(in),std::end(in),std::begin(m_ss));
    setState(State::ST_To_T); // (m_ss->m_rr)
    worker_ST_To_T(0);
    setState(State::Wait);
    std::copy(std::begin(m_rr),std::end(m_rr),std::begin(out));
}                                  

void HankelFFTRRVT::SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    // Transpose from SRST orientation to STSR orientation
    transpose(in.data(),m_ss.data(),m_nr,m_nst);
    // Compute transform over R coordindate
    setState(State::CMP_SR_To_R); // (m_ss->m_rs)
    workerCMP_SR_To_R(0);
    setState(State::Wait);
    // Transpose from STR orientation to RST orientation
    transpose(m_rs.data(),out.data(),m_nst,m_nr);
}

// Added
void HankelFFTRRVT::RST_To_SRST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out){
    // Transpose from RST orientation to STR orientation
    transpose(in.data(),m_rs.data(),m_nr,m_nst);
    // Compute transform over R coordinate
    setState(State::CMP_R_To_SR); // (m_rs->m_ss)
    workerCMP_R_To_SR(0);
    setState(State::Wait);
    // Transpose from STSR orientation to SRST orientation
    transpose(m_ss.data(),out.data(),m_nst,m_nr);
}


void HankelFFTRRVT::SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out)
{
    std::copy(std::begin(in),std::end(in),std::begin(m_ss));
    setState(State::ST_To_T); // (m_ss->m_rr)
    worker_ST_To_T(0);
    setState(State::Wait);
    std::copy(std::begin(m_rr),std::end(m_rr),std::begin(out));
}

void HankelFFTRRVT::SRT_To_SRST(const std::vector<double>& in, std::vector<dcmplx>& out)
{
    std::copy(std::begin(in),std::end(in),std::begin(m_rr));
    setState(State::T_To_ST); // (m_rr->m_ss)
    worker_T_To_ST(0);
    setState(State::Wait);
    std::copy(std::begin(m_ss),std::end(m_ss),std::begin(out));
}


void HankelFFTRRVT::RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out) 
{
    std::copy(std::begin(in),std::end(in),std::begin(m_rr));
    setState(State::T_To_ST); //(m_rr->m_ss)
    worker_T_To_ST(0);
    setState(State::Wait);
    std::copy(std::begin(m_ss),std::end(m_ss),std::begin(out));
}

void HankelFFTRRVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) 
{
    transpose(in.data(),m_ss.data(),m_nr,m_nst);
    setState(State::CMP_SR_To_R); // (m_ss->m_rs)
    workerCMP_SR_To_R(0);
    setState(State::Wait);

    transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
    setState(State::ST_To_T); // (m_ss->m_rr)
    worker_ST_To_T(0);
    setState(State::Wait);
    std::copy(std::begin(m_rr),std::end(m_rr),std::begin(out));
}

void HankelFFTRRVT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) 
{
    std::copy(std::begin(in),std::end(in),std::begin(m_rr));
    setState(State::T_To_ST); // (m_rr->m_ss)
    worker_T_To_ST(0);
    setState(State::Wait);

    transpose(m_ss.data(),m_rs.data(),m_nr,m_nst);

    setState(State::CMP_R_To_SR); // (m_rs->m_ss)
    workerCMP_R_To_SR(0);
    setState(State::Wait);

    transpose(m_ss.data(),out.data(),m_nst,m_nr);
}                                  

void HankelFFTRRVT::SRST_To_RCVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
{
    transpose(in.data(),m_ss.data(),m_nr,m_nst);
    setState(State::CMP_SR_To_R); // (m_ss->m_rs)
    workerCMP_SR_To_R(0);
    setState(State::Wait);

    transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
    setState(State::ST_To_CVT); // (m_ss->m_rs)
    worker_ST_To_CVT(0);
    setState(State::Wait);
    std::copy(std::begin(m_rs),std::end(m_rs),std::begin(out));
}

void HankelFFTRRVT::RCVT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
{
    std::copy(std::begin(in),std::end(in),std::begin(m_rs));
    setState(State::CVT_To_ST); // (m_rs->m_ss)
    worker_CVT_To_ST(0);
    setState(State::Wait);

    transpose(m_ss.data(),m_rs.data(),m_nr,m_nst);

    setState(State::CMP_R_To_SR); // (m_rs->m_ss)
    workerCMP_R_To_SR(0);
    setState(State::Wait);
    transpose(m_ss.data(),out.data(),m_nst,m_nr);
}                                  

void HankelFFTRRVT::RCVT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
{
    std::copy(std::begin(in),std::end(in),std::begin(m_rs));
    setState(State::CVT_To_ST); // (m_rs->m_ss)
    worker_CVT_To_ST(0);
    setState(State::Wait);
    std::copy(std::begin(m_ss),std::end(m_ss),std::begin(out));
}                              

void HankelFFTRRVT::RST_To_RCVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    std::copy(std::begin(in),std::end(in),std::begin(m_ss));
    setState(State::ST_To_CVT); // (m_ss->m_rs)
    worker_ST_To_CVT(0);
    setState(State::Wait);
    std::copy(std::begin(m_rs),std::end(m_rs),std::begin(out));
}                                  

void HankelFFTRRVT::worker_thread(unsigned tid)
{
  bool threads_done = false;
  while(!threads_done)
  {
      std::unique_lock lock{m_mut};
      m_cv.wait(lock,[this, tid]{return m_ready[tid-1];} );
      lock.unlock();

      if(m_state == State::T_To_ST) { // (m_rr->m_ss)
        worker_T_To_ST(tid);
        worker_wait(tid);
      }
      else if(m_state == State::ST_To_T){ // (m_ss->m_rr)
        worker_ST_To_T(tid);
        worker_wait(tid);
      }
      else if(m_state == State::CMP_R_To_SR){ // (m_rs->m_ss)
        workerCMP_R_To_SR(tid);
        worker_wait(tid);
      }
      else if(m_state == State::CMP_SR_To_R){ // (m_ss->m_rs)
        workerCMP_SR_To_R(tid);
        worker_wait(tid);
      }
      else if(m_state == State::R_To_SR){ // (m_rr->m_sr)
        worker_R_To_SR(tid);
        worker_wait(tid);
      }
      else if(m_state == State::SR_To_R){ // (m_sr->m_rr)
        worker_SR_To_R(tid);
        worker_wait(tid);
      }
      else if(m_state == State::CVT_To_ST){ // (m_rs->m_ss)
        worker_CVT_To_ST(tid);
        worker_wait(tid);
      }
      else if(m_state == State::ST_To_CVT){ // (m_ss->m_rs)
        worker_ST_To_CVT(tid);
        worker_wait(tid);
      }
      else if(m_state == State::Wait)
        continue;
      else if(m_state == State::Done)
        threads_done = true;
  }
}

void HankelFFTRRVT::worker_wait(unsigned tid)
{
    std::scoped_lock lock{m_mut};
    m_ready[tid - 1] = false;
    m_threads_processed++;
    if(m_threads_processed == m_thread.size()){
        m_threads_processed = 0;
        m_state = State::Wait;
        m_processed = true;
        m_cv.notify_all();
    }
}

void HankelFFTRRVT::setState(State state) {
    if(m_thread.empty())
        return;
    m_state = state;
    if(state == State::Wait){
        // Wait for threads to finish processing
        std::unique_lock lock{m_mut};
        m_cv.wait(lock,[this]{return m_processed;} );
        lock.unlock();
    } else if(state == State::Done){
        m_processed = true;
        for(auto item : m_ready)
            item = true;
        m_cv.notify_all();
    }
    else{
        // Activate state for thread processing
        m_processed = false;
        for(auto item : m_ready)
            item = true;
        // Tell threads to wake up and process, they are ready (m_ready)
        m_cv.notify_all();
    }
}

}