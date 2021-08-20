
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
        m_uRT(gridR.getNr()*gridT.getNt()),
        m_uRST(gridR.getNr()*gridT.getNst()),
        m_uSRST(gridR.getNr()*gridT.getNst())
    {
        m_nThreads = threads;
        m_nr = gridR.getNr();
        m_nt = gridT.getNt();
        m_nst = gridT.getNst();

        m_STATE = 0; // Do nothing
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
      ReadySTATE(5);
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

  void HankelPeriodicTransformRT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) 
  {
    std::copy(std::begin(in),std::end(in),std::begin(m_uRT));
    ReadySTATE(1);
    worker_RT_To_RST(0);
    ProcessedSTATE(0);
    // Transpose array
    for (unsigned int i = 0; i < m_nr; i++) 
        for (unsigned int j = 0; j < m_nst; j++)  
            m_uRST[m_nst*i + j] = m_uSRST[m_nr*j+i];

    ReadySTATE(2);
    worker_RST_To_SRST(0);
    ProcessedSTATE(0);
    // Transpose again
    for (unsigned int i = 0; i < m_nst; i++)  
        for (unsigned int j = 0; j < m_nr; j++) 
            out[m_nr*i + j] = m_uSRST[m_nr*j+i];
  }                                  

  void HankelPeriodicTransformRT::worker_RT_To_RST(int tid)
  {
      for (unsigned int i = tid; i < m_nr; i+=m_nThreads) 
          m_transformT[tid]->T_To_ST(m_uRT.data()+m_nt*i,m_uSRST.data()+m_nst*i);
  }

  void HankelPeriodicTransformRT::worker_RST_To_SRST(int tid)
  {
      for(unsigned int i = tid; i < m_nst; i+=m_nThreads)
          m_transformR[tid]->R_To_SR(m_uRST.data()+m_nr*i,m_uSRST.data()+m_nr*i);
  }

  void HankelPeriodicTransformRT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) 
  {
      // Transpose array
      for (unsigned int i = 0; i < m_nr; i++) 
          for (unsigned int j = 0; j < m_nst; j++)  
              m_uSRST[m_nst*i + j] = in[m_nr*j+i];

      ReadySTATE(3);
      worker_SRST_To_RST(0);
      ProcessedSTATE(0);

      // Transpose again 
      for(int i = 0; i < m_nr; i++)
          for(int j = 0; j < m_nst; j++)
              m_uSRST[m_nst*i+j] = m_uRST[m_nr*j+i];

      ReadySTATE(4);
      worker_RST_To_RT(0);
      ProcessedSTATE(0);
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

      if(m_STATE == 1) {
        worker_RT_To_RST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == 2){
        worker_RST_To_SRST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == 3){ 
        worker_SRST_To_RST(tid);
        worker_wait(tid);
      }
      else if(m_STATE == 4){
        worker_RST_To_RT(tid);
        worker_wait(tid);
      }
      else if(m_STATE == 0)
        continue;
      else if(m_STATE == 5)
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
          m_STATE = 0;
          m_processed = true;
          m_cv.notify_all();
      }
      m_mut.unlock();
  }


  void HankelPeriodicTransformRT::ReadySTATE(int beg_state)  
  {
    if(m_ready.size() > 0){
      m_STATE = beg_state;
      m_processed = false;
      {
        std::lock_guard<std::mutex> lock(m_mut);
        for(unsigned int i = 0; i < m_ready.size(); i++)
          m_ready[i] = true;
      }
      m_cv.notify_all();
    }
  }

  void HankelPeriodicTransformRT::ProcessedSTATE(int end_state)  
  {
    if(m_ready.size() > 0)
    {
      std::unique_lock<std::mutex> lock(m_mut);
      m_cv.wait(lock,[&]{return m_processed;} );
      m_STATE = end_state;
    }
  }






}
















