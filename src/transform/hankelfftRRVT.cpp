
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
        m_rs(m_nr*m_nst),
        m_sr(m_nr*m_nt),
        m_ss(m_nr*m_nst),
        m_thread(),
        m_ready()
    {
        m_STATE = State::Wait; // Do nothing
        m_THCOUNT = 0; 
        m_processed = true;

        for (unsigned m = 1; m < m_nThreads; m++)
        {
            m_ready.push_back(false);
            m_thread.push_back(std::thread(&HankelFFTRRVT::worker_thread,this,m));
        }

        for(unsigned i = 0; i < m_nThreads; i++) {
            FFTRVT* t_transform = new FFTRVT(gridT);
            HankelTransformR* r_transform = new HankelTransformR(gridR);
            m_transformT.push_back(t_transform);
            m_transformR.push_back(r_transform);
        }
  }

  HankelFFTRRVT::~HankelFFTRRVT(){
      ReadySTATE(State::Done);
      for (unsigned m = 0; m < m_thread.size(); m++)
          m_thread[m].join();
      m_thread.clear();
      for(auto item : m_transformT)
          delete item;
      for(auto item : m_transformR)
          delete item;
      m_transformT.clear();
      m_transformR.clear();
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
      ReadySTATE(State::R_To_SR); // (m_rr->m_sr)
      worker_R_To_SR(0);
      ProcessedSTATE(State::Wait);
      transpose(m_sr.data(),out.data(),m_nt,m_nr);
  }                                  

  void HankelFFTRRVT::SRT_To_RT(const std::vector<double>& in,std::vector<double>& out) 
  {
      transpose(in.data(),m_sr.data(),m_nr,m_nt);
      ReadySTATE(State::SR_To_R); // (m_sr->m_rr)
      worker_SR_To_R(0);
      ProcessedSTATE(State::Wait);
      transpose(m_rr.data(),out.data(),m_nt,m_nr);
  }


  void HankelFFTRRVT::RST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out)
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_ss));
      ReadySTATE(State::ST_To_T); // (m_ss->m_rr)
      worker_ST_To_T(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_rr),std::end(m_rr),std::begin(out));
  }                                  

  void HankelFFTRRVT::SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
  {
      // Transpose from SRST orientation to STSR orientation
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      // Compute transform over R coordindate
      ReadySTATE(State::CMP_SR_To_R); // (m_ss->m_rs)
      workerCMP_SR_To_R(0);
      ProcessedSTATE(State::Wait);
      // Transpose from STR orientation to RST orientation
      transpose(m_rs.data(),out.data(),m_nst,m_nr);
  }

  // Added
  void HankelFFTRRVT::RST_To_SRST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out){
      // Transpose from RST orientation to STR orientation
      transpose(in.data(),m_rs.data(),m_nr,m_nst);
      // Compute transform over R coordinate
      ReadySTATE(State::CMP_R_To_SR); // (m_rs->m_ss)
      workerCMP_R_To_SR(0);
      ProcessedSTATE(State::Wait);
      // Transpose from STSR orientation to SRST orientation
      transpose(m_ss.data(),out.data(),m_nst,m_nr);
  }


  void HankelFFTRRVT::SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out)
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_ss));
      ReadySTATE(State::ST_To_T); // (m_ss->m_rr)
      worker_ST_To_T(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_rr),std::end(m_rr),std::begin(out));
  }

  void HankelFFTRRVT::SRT_To_SRST(const std::vector<double>& in, std::vector<dcmplx>& out)
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_rr));
      ReadySTATE(State::T_To_ST); // (m_rr->m_ss)
      worker_T_To_ST(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_ss),std::end(m_ss),std::begin(out));
  }


  void HankelFFTRRVT::RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_rr));
      ReadySTATE(State::T_To_ST); //(m_rr->m_ss)
      worker_T_To_ST(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_ss),std::end(m_ss),std::begin(out));
  }

  void HankelFFTRRVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) 
  {
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      ReadySTATE(State::CMP_SR_To_R); // (m_ss->m_rs)
      workerCMP_SR_To_R(0);
      ProcessedSTATE(State::Wait);

      transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
      ReadySTATE(State::ST_To_T); // (m_ss->m_rr)
      worker_ST_To_T(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_rr),std::end(m_rr),std::begin(out));
  }

  void HankelFFTRRVT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_rr));
      ReadySTATE(State::T_To_ST); // (m_rr->m_ss)
      worker_T_To_ST(0);
      ProcessedSTATE(State::Wait);

      transpose(m_ss.data(),m_rs.data(),m_nr,m_nst);

      ReadySTATE(State::CMP_R_To_SR); // (m_rs->m_ss)
      workerCMP_R_To_SR(0);
      ProcessedSTATE(State::Wait);

      transpose(m_ss.data(),out.data(),m_nst,m_nr);
  }                                  

  void HankelFFTRRVT::SRST_To_RCVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      ReadySTATE(State::CMP_SR_To_R); // (m_ss->m_rs)
      workerCMP_SR_To_R(0);
      ProcessedSTATE(State::Wait);

      transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
      ReadySTATE(State::ST_To_CVT); // (m_ss->m_rs)
      worker_ST_To_CVT(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_rs),std::end(m_rs),std::begin(out));
  }

  void HankelFFTRRVT::RCVT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_rs));
      ReadySTATE(State::CVT_To_ST); // (m_rs->m_ss)
      worker_CVT_To_ST(0);
      ProcessedSTATE(State::Wait);

      transpose(m_ss.data(),m_rs.data(),m_nr,m_nst);

      ReadySTATE(State::CMP_R_To_SR); // (m_rs->m_ss)
      workerCMP_R_To_SR(0);
      ProcessedSTATE(State::Wait);
      transpose(m_ss.data(),out.data(),m_nst,m_nr);
  }                                  

  void HankelFFTRRVT::RCVT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_rs));
      ReadySTATE(State::CVT_To_ST); // (m_rs->m_ss)
      worker_CVT_To_ST(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_ss),std::end(m_ss),std::begin(out));
  }                              

  void HankelFFTRRVT::RST_To_RCVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_ss));
      ReadySTATE(State::ST_To_CVT); // (m_ss->m_rs)
      worker_ST_To_CVT(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_rs),std::end(m_rs),std::begin(out));
  }                                  

  void HankelFFTRRVT::worker_thread(unsigned tid)
  {
    bool threads_done = false;
    while(!threads_done)
    {
      {
        std::unique_lock<std::mutex> lock(m_mut);
        m_cv.wait(lock,[&]{return m_ready[tid-1];} );
      }

      if(m_STATE == State::T_To_ST) { // (m_rr->m_ss)
        worker_T_To_ST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::ST_To_T){ // (m_ss->m_rr)
        worker_ST_To_T(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::CMP_R_To_SR){ // (m_rs->m_ss)
        workerCMP_R_To_SR(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::CMP_SR_To_R){ // (m_ss->m_rs)
        workerCMP_SR_To_R(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::R_To_SR){ // (m_rr->m_sr)
        worker_R_To_SR(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::SR_To_R){ // (m_sr->m_rr)
        worker_SR_To_R(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::CVT_To_ST){ // (m_rs->m_ss)
        worker_CVT_To_ST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::ST_To_CVT){ // (m_ss->m_rs)
        worker_ST_To_CVT(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::Wait)
        continue;
      else if(m_STATE == State::Done)
        threads_done = true;
    }
  }

  void HankelFFTRRVT::worker_wait(unsigned tid)
  {
      m_mut.lock();
      m_ready[tid - 1] = false;
      m_THCOUNT++;
      if(m_THCOUNT == m_nThreads - 1){
          m_THCOUNT = 0;
          m_STATE = State::Wait;
          m_processed = true;
          m_cv.notify_all();
      }
      m_mut.unlock();
  }


  void HankelFFTRRVT::ReadySTATE(State state)  
  {
      if(m_nThreads < 2)
          return;

      m_STATE = state;
      m_processed = false;
      {
          std::lock_guard<std::mutex> lock(m_mut);
          for(unsigned i = 0; i < m_ready.size(); i++)
              m_ready[i] = true;
      }
      m_cv.notify_all();
  }

  void HankelFFTRRVT::ProcessedSTATE(State state)  
  {
      if(m_nThreads < 2)
          return;
      std::unique_lock<std::mutex> lock(m_mut);
      m_cv.wait(lock,[&]{return m_processed;} );
      m_STATE = state;
  }










/*
  HankelFFTRRVT::HankelFFTRRVT(const BesselRootGridR& gridR,\
          const UniformGridRVT& gridT,unsigned threads) :
        m_nr(gridR.getNr()),
        m_nt(gridT.getNt()),
        m_nst(gridT.getNst()),
        m_threads(threads),
        m_rr(m_nr*m_nt),
        m_rs(m_nr*m_nst),
        m_sr(m_nr*m_nt),
        m_ss(m_nr*m_nst),
        m_transformT(threads,nullptr),
        m_transformR(threads,nullptr)
    {
        for(auto i = 0; i < threads; i++){
            m_transformT[i] = new FFTRVT(gridT);
            m_transformR[i] = new HankelTransformR(gridR);
        }
    }

  HankelFFTRRVT::~HankelFFTRRVT(){
      for(auto item : m_transformT)
          delete item;
      for(auto item : m_transformR)
          delete item;
  }

  void HankelFFTRRVT::worker_T_To_ST(unsigned tid,const double* in,dcmplx* out){
      for(unsigned i = tid; i < m_nr; i+= m_threads)
          m_transformT[tid]->T_To_ST(in+m_nt*i,out+m_nst*i);
  }

  void HankelFFTRRVT::worker_ST_To_T(unsigned tid,const dcmplx* in,double* out)
  {
      for (unsigned i = tid; i < m_nr; i+= m_threads)
          m_transformT[tid]->ST_To_T(in+m_nst*i,out+m_nt*i);
  }

  void HankelFFTRRVT::worker_CVT_To_ST(unsigned tid,const dcmplx* in,dcmplx* out){
      for(unsigned i = tid; i < m_nr; i+= m_threads)
          m_transformT[tid]->CVT_To_ST(in+m_nt*i,out+m_nst*i);
  }

  void HankelFFTRRVT::worker_ST_To_CVT(unsigned tid,const dcmplx* in,dcmplx* out)
  {
      for (unsigned i = tid; i < m_nr; i+= m_threads)
          m_transformT[tid]->ST_To_CVT(in+m_nst*i,out+m_nt*i);
  }

  void HankelFFTRRVT::worker_R_To_SR(unsigned tid,const double* in,double* out){
      for(unsigned i = tid; i < m_nt; i+= m_threads)
          m_transformR[tid]->R_To_SR(in+m_nr*i,out+m_nr*i);
  }

  void HankelFFTRRVT::worker_SR_To_R(unsigned tid,const double* in,double* out)
  {
      for (unsigned j = tid; j < m_nt; j+=m_threads)  
          m_transformR[tid]->SR_To_R(in+m_nr*j,out+m_nr*j);
  }

  void HankelFFTRRVT::workerCMP_R_To_SR(unsigned tid,const dcmplx* in,dcmplx* out){
      for(unsigned i = tid; i < m_nst; i+= m_threads)
          m_transformR[tid]->R_To_SR(in+m_nr*i,out+m_nr*i);
  }

  void HankelFFTRRVT::workerCMP_SR_To_R(unsigned tid,const dcmplx* in, dcmplx* out)
  {
      for (unsigned j = tid; j < m_nst; j+=m_threads)  
          m_transformR[tid]->SR_To_R(in+m_nr*j,out+m_nr*j);
  }

  void HankelFFTRRVT::wait_for_workers(std::vector<std::thread>& workers)
  {
      for(auto& worker : workers)
          worker.join();
      workers.clear();
  }

  void HankelFFTRRVT::RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_T_To_ST,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  

  void HankelFFTRRVT::CVRT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_CVT_To_ST,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                              


  void HankelFFTRRVT::RST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_ST_To_T,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  

  void HankelFFTRRVT::RST_To_CVRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_ST_To_CVT,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  


  void HankelFFTRRVT::RST_To_SRST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out){
      // Transpose from RST orientation to STR orientation
      transpose(in.data(),m_rs.data(),m_nr,m_nst);
      // Compute transform over R coordinate
      HankelFFTRRVT::hSTR_To_STSR(m_rs,m_ss);
      // Transpose from STSR orientation to SRST orientation
      transpose(m_ss.data(),out.data(),m_nst,m_nr);
  }

  void HankelFFTRRVT::SRST_To_RST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out){
      // Transpose from SRST orientation to STSR orientation
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      // Compute transform over R coordinate
      HankelFFTRRVT::hSTSR_To_STR(m_ss,m_rs);
      // Transpose from STR orientation to RST orientation
      transpose(m_rs.data(),out.data(),m_nst,m_nr);
  }

  void HankelFFTRRVT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      HankelFFTRRVT::RT_To_RST(in,m_rs);
      transpose(m_rs.data(),m_ss.data(),m_nr,m_nst);
      HankelFFTRRVT::hSTR_To_STSR(m_ss,m_rs);
      transpose(m_rs.data(),out.data(),m_nst,m_nr);
  }                                  

  void HankelFFTRRVT::CVRT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      HankelFFTRRVT::CVRT_To_RST(in,m_rs);
      transpose(m_rs.data(),m_ss.data(),m_nr,m_nst);
      HankelFFTRRVT::hSTR_To_STSR(m_ss,m_rs);
      transpose(m_rs.data(),out.data(),m_nst,m_nr);
  }                                  


  void HankelFFTRRVT::RT_To_SRT(const std::vector<double>& in,std::vector<double>& out) 
  {
      transpose(in.data(),m_rr.data(),m_nr,m_nt);
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_R_To_SR,this,\
                      tid,m_rr.data(),m_sr.data()));
      }
      wait_for_workers(workers);
      transpose(m_sr.data(),out.data(),m_nt,m_nr);
  }

  void HankelFFTRRVT::SRT_To_RT(const std::vector<double>& in,std::vector<double>& out) 
  {
      transpose(in.data(),m_sr.data(),m_nr,m_nt);
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_SR_To_R,this,\
                      tid,m_sr.data(),m_rr.data()));
      }
      wait_for_workers(workers);
      transpose(m_rr.data(),out.data(),m_nt,m_nr);
  }

  void HankelFFTRRVT::SRT_To_SRST(const std::vector<double>& in, std::vector<dcmplx>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_T_To_ST,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }

  void HankelFFTRRVT::SRST_To_SRT(const std::vector<dcmplx>& in, std::vector<double>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::worker_ST_To_T,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }

  void HankelFFTRRVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) 
  {
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      HankelFFTRRVT::hSTSR_To_STR(m_ss,m_rs);
      transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
      HankelFFTRRVT::RST_To_RT(m_ss,out);
  }

  void HankelFFTRRVT::SRST_To_CVRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      HankelFFTRRVT::hSTSR_To_STR(m_ss,m_rs);
      transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
      HankelFFTRRVT::RST_To_CVRT(m_ss,out);
  }

  void HankelFFTRRVT::hSTR_To_STSR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::workerCMP_R_To_SR,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  

  void HankelFFTRRVT::hSTSR_To_STR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRRVT::workerCMP_SR_To_R,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }
  */



  /* 
   *
   *
 
  // Attempt at using a third-party thread pool:
  // -> some speed-up with multiple threads but slower in general and doesnt scale well

  HankelFFTRRVTc::HankelFFTRRVTc(const BesselRootGridR& gridR,const UniformGridT& gridT) :
        m_nr(gridR.getNr()),
        m_nt(gridT.getNt()),
        m_nst(gridT.getNst()),
        m_rr(m_nr*m_nt),
        m_rs(m_nr*m_nst),
        m_sr(m_nr*m_nt),
        m_ss(m_nr*m_nst),
        m_gridR(new BesselRootGridR(gridR)),
        m_gridT(new UniformGridT(gridT)),
        m_fftblt(new FFTRVT(gridT)),
        m_hankel(new HankelTransformR(gridR)),
        m_transformT(),
        m_transformR(),
        m_pool(nullptr) { } 

  HankelFFTRRVTc::~HankelFFTRRVTc(){


      delete m_gridR;
      delete m_gridT;
      delete m_fftblt;
      delete m_hankel;
      for(auto item : m_transformT)
          delete item.second;
      for(auto item : m_transformR)
          delete item.second;
      m_transformT.clear();
      m_transformR.clear();
  }

  bool HankelFFTRRVTc::setThreadPool(thread_pool* pool)
  {
      if(!pool)
          return false;

      m_pool = pool;
      std::vector<std::thread::id> ids = m_pool->get_thread_ids();
      for(auto i = 0; i < m_pool->get_thread_count(); i++){
          FFTRVT* t_transform = new FFTRVT(*m_gridT);
          HankelTransformR* r_transform = new HankelTransformR(*m_gridR);
          m_transformT.insert(std::pair(ids[i],t_transform));
          m_transformR.insert(std::pair(ids[i],r_transform));
      }
      return true;
  }

  void worker_RT_To_RST(const double* in,dcmplx* out,std::map<std::thread::id,FFTRVT*> fftblt_map)
  {
      std::thread::id this_id = std::this_thread::get_id();
      fftblt_map[this_id]->T_To_ST(in,out);
  }

  void worker_STR_To_STSR(const dcmplx* in,dcmplx* out,std::map<std::thread::id,\
          HankelTransformR*> hankel_map)
  {
      std::thread::id this_id = std::this_thread::get_id();
      hankel_map[this_id]->R_To_SR(in,out);
  }

  void HankelFFTRRVTc::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      if(!m_pool){
          for(unsigned i = 0; i < m_nr; i++)
              m_fftblt->T_To_ST(in.data()+m_nt*i,m_ss.data()+m_nst*i);
          // transpose
          for(unsigned i = 0; i < m_nr; i++)
              for (unsigned j = 0; j < m_nst; j++){
                  m_rs[j*m_nr+i] = m_ss[i*m_nst+j];
          }
          for(unsigned i = 0; i < m_nst; i++)
              m_hankel->R_To_SR(m_rs.data()+m_nr*i,m_ss.data()+m_nr*i);
                // transpose again
          for(unsigned i = 0; i < m_nr; i++)
              for(unsigned j = 0; j < m_nst; j++)
                  out[i*m_nst+j] = m_ss[m_nr*j + i]; 
      } else {
          for(unsigned i = 0; i < m_nr; i++){
              m_pool->push_task(std::bind(spida::worker_RT_To_RST,\
                      in.data()+m_nt*i,m_ss.data()+m_nst*i,m_transformT));
          }
          m_pool->wait_for_tasks();

          auto transpose = [] (const uint32_t& iS,const uint32_t& iF,\
                  uint32_t jS,uint32_t jF,dcmplx* A,dcmplx* B) {
              for(uint32_t i = iS; i < iF; i++)
                  for (uint32_t j = jS; j < jF; j++)
                      B[j*iF+i] = A[i*jF+j];
          };

          m_pool->parallelize_loop(0,m_nr,std::bind(transpose,std::placeholders::_1,
                      std::placeholders::_2,0,m_nst,m_ss.data(),m_rs.data()));

          m_pool->wait_for_tasks();

          for(unsigned i = 0; i < m_nst; i++){
              m_pool->push_task(std::bind(spida::worker_STR_To_STSR,\
                      m_rs.data()+m_nr*i,m_ss.data()+m_nr*i,m_transformR));
          }
          m_pool->wait_for_tasks();


          m_pool->parallelize_loop(0,m_nr,std::bind(transpose,std::placeholders::_1,
                      std::placeholders::_2,0,m_nst,m_ss.data(),out.data()));
          m_pool->wait_for_tasks();


      }
  }                                  

  void HankelFFTRRVTc::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) 
  {
      for(unsigned i = 0; i < m_nr; i++)
            for(unsigned j = 0; j < m_nst; j++)
                m_ss[m_nr*j+i] = in[m_nst*i + j]; 

      if(!m_pool){
          for (unsigned j = 0; j < m_nst; j++)  
              m_hankel->SR_To_R(m_ss.data()+m_nr*j,m_rs.data()+m_nr*j);
      }

      for(unsigned i = 0; i < m_nr; i++)
          for(unsigned j = 0; j < m_nst; j++)
              m_ss[m_nst*i+j] = m_rs[m_nr*j + i]; 

      if(!m_pool){
          for (unsigned j = 0; j < m_nst; j++)  
              m_fftblt->ST_To_T(m_ss.data()+m_nst*j,out.data()+m_nt*j);
      }
  }

  */




}
















