
#include <algorithm>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "spida/transform/hankelperiodicRT.h"
#include "spida/transform/periodicT.h"
#include "spida/transform/hankelR.h"
#include "spida/grid/besselR.h" 
#include "spida/grid/uniformT.h" 

namespace spida {

    HankelPeriodicTransformRT::HankelPeriodicTransformRT(const BesselRootGridR& gridR,const UniformGridT& gridT,int threads) :
        m_nr(gridR.getNr()),
        m_nt(gridT.getNt()),
        m_nst(gridT.getNst()),
        m_nThreads(threads),
        m_uRT(m_nr*m_nt),
        m_uRST(m_nr*m_nst),
        m_uSRT(m_nr*m_nt),
        m_uSRST(m_nr*m_nst),
        m_thread(),
        m_ready()
    {
        m_STATE = State::Wait; // Do nothing
        m_THCOUNT = 0; 
        m_processed = true;

        for (int m = 1; m < m_nThreads; m++)
        {
            m_ready.push_back(false);
            m_thread.push_back(std::thread(&HankelPeriodicTransformRT::worker_thread,this,m));
        }

        for(int i = 0; i < m_nThreads; i++) {
            PeriodicTransformT* t_transform = new PeriodicTransformT(gridT);
            HankelTransformR* r_transform = new HankelTransformR(gridR);
            m_transformT.push_back(t_transform);
            m_transformR.push_back(r_transform);
        }
  }

  HankelPeriodicTransformRT::~HankelPeriodicTransformRT(){
      ReadySTATE(State::Done);
      for (unsigned int m = 0; m < m_thread.size(); m++)
          m_thread[m].join();
      m_thread.clear();
      for(auto item : m_transformT)
          delete item;
      for(auto item : m_transformR)
          delete item;
      m_transformT.clear();
      m_transformR.clear();
  }

  void HankelPeriodicTransformRT::RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_uRT));
      ReadySTATE(State::RT_To_RST);
      worker_RT_To_RST(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_uSRST),std::end(m_uSRST),std::begin(out));
  }                                  

  void HankelPeriodicTransformRT::RT_To_SRT(const std::vector<double>& in,std::vector<double>& out) 
  {
//      for(unsigned int i = 0; i < m_nr; i++)
//          for (unsigned int j = 0; j < m_nt; j++){
//              m_uRT[j*m_nr+i] = in[i*m_nt+j];
//      }
//      for(unsigned int j = 0; j < m_nt; j++)
//          m_transformR[0]->R_To_SR(m_uRT.data()+j*m_nr,m_uSRT.data()+j*m_nr);
//
//      for(unsigned int i = 0; i < m_nr; i++)
//          for(unsigned int j = 0; j < m_nt; j++)
//              out[i*m_nt+j] = m_uSRT[m_nr*j + i]; 

      for(unsigned int i = 0; i < m_nr; i++)
          for (unsigned int j = 0; j < m_nt; j++){
              m_uRT[j*m_nr+i] = in[i*m_nt+j];
      }
      ReadySTATE(State::RT_To_SRT);
      worker_RT_To_SRT(0);
      ProcessedSTATE(State::Wait);
      for(unsigned int i = 0; i < m_nr; i++)
          for(unsigned int j = 0; j < m_nt; j++)
              out[i*m_nt+j] = m_uSRT[m_nr*j + i]; 
  }                                  

  void HankelPeriodicTransformRT::SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
  {
      return;
  }

  void HankelPeriodicTransformRT::SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out)
  {
      return;
  }

  void HankelPeriodicTransformRT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
      std::copy(std::begin(in),std::end(in),std::begin(m_uRT));
      // Transform from T to ST
      ReadySTATE(State::RT_To_RST);
      worker_RT_To_RST(0);
      ProcessedSTATE(State::Wait);

      // Transpose array -> R data ready for transform
      for (unsigned int j = 0; j < m_nst; j++)  
          for (unsigned int i = 0; i < m_nr; i++) 
              m_uRST[m_nr*i + j] = m_uSRST[i*m_nst + j];
  
      // Transform from R to SR
      ReadySTATE(State::RST_To_SRST);
      worker_RST_To_SRST(0);
      ProcessedSTATE(State::Wait);
      // Transpose again -> out is in SR-ST row-column form
      for (unsigned int i = 0; i < m_nst; i++)  
          for (unsigned int j = 0; j < m_nr; j++) 
              out[m_nr*i + j] = m_uSRST[m_nst*j+i];
  }                                  

  void HankelPeriodicTransformRT::worker_RT_To_RST(int tid)
  {
      for (unsigned int i = tid; i < m_nr; i+=m_nThreads) 
          m_transformT[tid]->T_To_ST(m_uRT.data()+m_nt*i,m_uSRST.data()+m_nst*i);
  }

  void HankelPeriodicTransformRT::worker_RST_To_SRST(int tid)
  {
      for(unsigned int i = tid; i < m_nst; i+=m_nThreads)
          m_transformR[tid]->R_To_SR(reinterpret_cast<const dcmplx*>(m_uRST.data()+m_nr*i),\
                  m_uSRST.data()+m_nr*i);
  }

  void HankelPeriodicTransformRT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) 
  {
      // Transpose array
      for (unsigned int i = 0; i < m_nr; i++) 
          for (unsigned int j = 0; j < m_nst; j++)  
              m_uSRST[m_nst*i + j] = in[m_nr*j+i];

      ReadySTATE(State::SRST_To_RST);
      worker_SRST_To_RST(0);
      ProcessedSTATE(State::Wait);
      // Transpose again 
      for(int i = 0; i < m_nr; i++)
          for(int j = 0; j < m_nst; j++)
              m_uSRST[m_nst*i+j] = m_uRST[m_nr*j+i];

      ReadySTATE(State::RST_To_RT);
      worker_RST_To_RT(0);
      ProcessedSTATE(State::Wait);
      std::copy(std::begin(m_uRT),std::end(m_uRT),std::begin(out));
  }

  void HankelPeriodicTransformRT::worker_SRST_To_RST(int tid)
  {
      for (unsigned int i = tid; i < m_nst; i+=m_nThreads)  
          m_transformR[tid]->SR_To_R(m_uSRST.data()+m_nr*i,m_uRST.data()+m_nr*i);
  }

  void HankelPeriodicTransformRT::worker_RST_To_RT(int tid)
  {
      for (unsigned int i = tid; i < m_nr; i+= m_nThreads)
          m_transformT[tid]->ST_To_T(m_uSRST.data()+m_nst*i,m_uRT.data()+m_nst*i);
  }


  void HankelPeriodicTransformRT::worker_thread(int tid)
  {
    bool threads_done = false;
    while(!threads_done)
    {
      {
        std::unique_lock<std::mutex> lock(m_mut);
        m_cv.wait(lock,[&]{return m_ready[tid-1];} );
      }

      if(m_STATE == State::RT_To_RST) {
        worker_RT_To_RST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::RST_To_SRST){
        worker_RST_To_SRST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::SRST_To_RST){ 
        worker_SRST_To_RST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::RST_To_RT){
        worker_RST_To_RT(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::RT_To_SRT){
        worker_RT_To_SRT(tid);
        worker_wait(tid);
      }
      else if(m_STATE == State::SRST_To_SRT){
        worker_SRST_To_SRT(tid);
        worker_wait(tid);   
      }
      else if(m_STATE == State::Wait)
        continue;
      else if(m_STATE == State::Done)
        threads_done = true;
    }
  }

  void HankelPeriodicTransformRT::worker_wait(int tid)
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

  void HankelPeriodicTransformRT::worker_RT_To_SRT(int tid)
  {
      for(unsigned int j = tid; j < m_nt; j+=m_nThreads)
          m_transformR[tid]->R_To_SR(m_uRT.data()+j*m_nr,m_uSRT.data()+j*m_nr);
  }

  void HankelPeriodicTransformRT::worker_SRST_To_SRT(int tid)
  {
      // not implemented
      return;
  }


  void HankelPeriodicTransformRT::ReadySTATE(State state)  
  {
      if(m_nThreads < 2)
          return;

      m_STATE = state;
      m_processed = false;
      {
          std::lock_guard<std::mutex> lock(m_mut);
          for(unsigned int i = 0; i < m_ready.size(); i++)
              m_ready[i] = true;
      }
      m_cv.notify_all();
  }

  void HankelPeriodicTransformRT::ProcessedSTATE(State state)  
  {
      if(m_nThreads < 2)
          return;
      std::unique_lock<std::mutex> lock(m_mut);
      m_cv.wait(lock,[&]{return m_processed;} );
      m_STATE = state;
  }






}
















