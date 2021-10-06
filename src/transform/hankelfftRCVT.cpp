
#include <algorithm>
#include <vector>
#include <thread>
#include "spida/transform/hankelfftRCVT.h"
#include "spida/transform/fftCVT.h"
#include "spida/transform/hankelR.h"
#include "spida/grid/besselR.h" 
#include "spida/grid/uniformCVT.h" 
#include "spida/helper/funcs.hpp"

namespace spida {

  HankelFFTRCVT::HankelFFTRCVT(const BesselRootGridR& gridR,\
          const UniformGridCVT& gridT,unsigned threads) :
        m_nr(gridR.getNr()),
        m_nt(gridT.getNt()),
        m_nst(gridT.getNst()),
        m_threads(threads),
        m_rs(m_nr*m_nst),
        m_ss(m_nr*m_nst),
        m_rr(m_nr*m_nt),
        m_sr(m_nr*m_nt),
        m_transformT(threads,nullptr),
        m_transformR(threads,nullptr)
    {
        for(auto i = 0; i < threads; i++){
            m_transformT[i] = new FFTCVT(gridT);
            m_transformR[i] = new HankelTransformR(gridR);
        }
    }

  HankelFFTRCVT::~HankelFFTRCVT(){
      for(auto item : m_transformT)
          delete item;
      for(auto item : m_transformR)
          delete item;
  }

  void HankelFFTRCVT::worker_T_To_ST(unsigned tid,const dcmplx* in,dcmplx* out){
      for(unsigned i = tid; i < m_nr; i+= m_threads)
          m_transformT[tid]->T_To_ST(in+m_nt*i,out+m_nst*i);
  }

  void HankelFFTRCVT::worker_ST_To_T(unsigned tid,const dcmplx* in,dcmplx* out)
  {
      for (unsigned i = tid; i < m_nr; i+= m_threads)
          m_transformT[tid]->ST_To_T(in+m_nst*i,out+m_nt*i);
  }

  void HankelFFTRCVT::worker_R_To_SR(unsigned tid,const dcmplx* in,dcmplx* out){
      for(unsigned i = tid; i < m_nt; i+= m_threads)
          m_transformR[tid]->R_To_SR(in+m_nr*i,out+m_nr*i);
  }

  void HankelFFTRCVT::worker_SR_To_R(unsigned tid,const dcmplx* in,dcmplx* out)
  {
      for (unsigned j = tid; j < m_nt; j+=m_threads)  
          m_transformR[tid]->SR_To_R(in+m_nr*j,out+m_nr*j);
  }

  void HankelFFTRCVT::workerST_R_To_SR(unsigned tid,const dcmplx* in,dcmplx* out){
      for(unsigned i = tid; i < m_nst; i+= m_threads)
          m_transformR[tid]->R_To_SR(in+m_nr*i,out+m_nr*i);
  }

  void HankelFFTRCVT::workerST_SR_To_R(unsigned tid,const dcmplx* in, dcmplx* out)
  {
      for (unsigned j = tid; j < m_nst; j+=m_threads)  
          m_transformR[tid]->SR_To_R(in+m_nr*j,out+m_nr*j);
  }

  void HankelFFTRCVT::wait_for_workers(std::vector<std::thread>& workers)
  {
      for(auto& worker : workers)
          worker.join();
      workers.clear();
  }

  void HankelFFTRCVT::RT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_T_To_ST,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  

  void HankelFFTRCVT::RST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_ST_To_T,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  

  void HankelFFTRCVT::RST_To_SRST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out){
      // Transpose from RST orientation to STR orientation
      transpose(in.data(),m_rs.data(),m_nr,m_nst);
      // Compute transform over R coordinate
      HankelFFTRCVT::hSTR_To_STSR(m_rs,m_ss);
      // Transpose from STSR orientation to SRST orientation
      transpose(m_ss.data(),out.data(),m_nst,m_nr);
  }

  void HankelFFTRCVT::hSTR_To_STSR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_R_To_SR,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }                                  

  void HankelFFTRCVT::SRST_To_RST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out){
      // Transpose from SRST orientation to STSR orientation
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      // Compute transform over R coordinate
      HankelFFTRCVT::hSTSR_To_STR(m_ss,m_rs);
      // Transpose from STR orientation to RST orientation
      transpose(m_rs.data(),out.data(),m_nst,m_nr);
  }

  void HankelFFTRCVT::hSTSR_To_STR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_SR_To_R,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }

  void HankelFFTRCVT::RT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      HankelFFTRCVT::RT_To_RST(in,m_rs);
      transpose(m_rs.data(),m_ss.data(),m_nr,m_nst);
      HankelFFTRCVT::hSTR_To_STSR(m_ss,m_rs);
      transpose(m_rs.data(),out.data(),m_nst,m_nr);
  }                                  

  void HankelFFTRCVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      transpose(in.data(),m_ss.data(),m_nr,m_nst);
      HankelFFTRCVT::hSTSR_To_STR(m_ss,m_rs);
      transpose(m_rs.data(),m_ss.data(),m_nst,m_nr);
      HankelFFTRCVT::RST_To_RT(m_ss,out);
  }


  void HankelFFTRCVT::RT_To_SRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      transpose(in.data(),m_rr.data(),m_nr,m_nt);
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_R_To_SR,this,\
                      tid,m_rr.data(),m_sr.data()));
      }
      wait_for_workers(workers);
      transpose(m_sr.data(),out.data(),m_nt,m_nr);
  }

  void HankelFFTRCVT::SRT_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
  {
      transpose(in.data(),m_sr.data(),m_nr,m_nt);
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_SR_To_R,this,\
                      tid,m_sr.data(),m_rr.data()));
      }
      wait_for_workers(workers);
      transpose(m_rr.data(),out.data(),m_nt,m_nr);
  }

  void HankelFFTRCVT::SRT_To_SRST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_T_To_ST,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }

  void HankelFFTRCVT::SRST_To_SRT(const std::vector<dcmplx>& in, std::vector<dcmplx>& out)
  {
      std::vector<std::thread> workers;
      for(auto tid = 0; tid < m_threads; tid++){
          workers.push_back(std::thread(&HankelFFTRCVT::worker_ST_To_T,this,\
                      tid,in.data(),out.data()));
      }
      wait_for_workers(workers);
  }








}









