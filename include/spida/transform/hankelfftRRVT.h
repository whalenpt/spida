// hankelfftRRVT.h
// Class for 2D functions of radius (R) and band-limited time (RVT)
#pragma once

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "spida/helper/constants.h"

namespace spida{

class FFTRVT;
class HankelTransformR;
class BesselRootGridR;
class UniformGridRVT;


class HankelFFTRRVT 
{
    public:
        explicit HankelFFTRRVT(const BesselRootGridR& gridR,const UniformGridRVT& gridT,unsigned threads=1);
        ~HankelFFTRRVT();
        HankelFFTRRVT()=delete;
        HankelFFTRRVT(const HankelFFTRRVT& sp)=delete;
        HankelFFTRRVT& operator=(const HankelFFTRRVT& sp)=delete;

        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out);  // tested - success 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); // tested - success
        void RCVT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_RCVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);

        void RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out);  // tested - success
        void RST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out);  // tested - success

        void RT_To_SRT(const std::vector<double>& in,std::vector<double>& out);  // tested - success
        void SRT_To_RT(const std::vector<double>& in,std::vector<double>& out);  // tested - success
        void RCVT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void RST_To_RCVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 

        void RST_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); // tested - success
        void SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); // tested - success

        void SRT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); // tested - success
        void SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out); // tested - success

        enum class State{Wait,\
            T_To_ST,\
            ST_To_T,\
            CMP_R_To_SR,\
            CMP_SR_To_R,\
            R_To_SR,\
            SR_To_R,\
            CVT_To_ST,\
            ST_To_CVT,\
            Done};

    private:
        unsigned m_nr, m_nt, m_nst;
        unsigned m_nThreads;

        std::vector<double> m_rr, m_sr;
        std::vector<dcmplx> m_rs, m_ss;

        std::vector<std::thread> m_thread;
        std::vector<FFTRVT*> m_transformT;
        std::vector<HankelTransformR*> m_transformR;

        State m_STATE;
        unsigned m_THCOUNT;
        std::vector<bool> m_ready;
        bool m_processed;
        std::mutex m_mut;
        std::condition_variable m_cv;

        void ReadySTATE(State state);
        void ProcessedSTATE(State state);

        void worker_thread(unsigned id);
        void worker_wait(unsigned id);
        void worker_T_To_ST(unsigned tid);
        void worker_ST_To_T(unsigned tid);
        void worker_ST_To_CVT(unsigned tid);
        void worker_CVT_To_ST(unsigned tid);
        void worker_R_To_SR(unsigned tid);
        void worker_SR_To_R(unsigned tid);
        void workerCMP_R_To_SR(unsigned tid);
        void workerCMP_SR_To_R(unsigned tid);
};

/*
class HankelFFTRRVT 
{
    public:
        explicit HankelFFTRRVT(const BesselRootGridR& gridR,const UniformGridRVT& gridT,unsigned threads = 1);
        ~HankelFFTRRVT();
        HankelFFTRRVT()=delete;
        HankelFFTRRVT(const HankelFFTRRVT& sp)=delete;
        HankelFFTRRVT& operator=(const HankelFFTRRVT& sp)=delete;

        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 
        void CVRT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_CVRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);

        void RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void CVRT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void RST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 
        void RST_To_CVRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 

        void RST_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);

        void RT_To_SRT(const std::vector<double>& in,std::vector<double>& out); 
        void SRT_To_RT(const std::vector<double>& in,std::vector<double>& out); 

        void SRT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out);
        void SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out);

    private:
        int m_nr;
        int m_nt;
        int m_nst;
        int m_threads;

        std::vector<double> m_rr;
        std::vector<dcmplx> m_rs;
        std::vector<double> m_sr;
        std::vector<dcmplx> m_ss;

        std::vector<FFTRVT*> m_transformT;
        std::vector<HankelTransformR*> m_transformR;

        void worker_T_To_ST(unsigned tid,const double* in,dcmplx* out);
        void worker_ST_To_T(unsigned tid,const dcmplx* in,double* out);

        void worker_ST_To_CVT(unsigned tid,const dcmplx* in,dcmplx* out);
        void worker_CVT_To_ST(unsigned tid,const dcmplx* in,dcmplx* out);

        void worker_R_To_SR(unsigned tid,const double* in,double* out);
        void worker_SR_To_R(unsigned tid,const double* in,double* out);
        void workerCMP_R_To_SR(unsigned tid,const dcmplx* in,dcmplx* out);
        void workerCMP_SR_To_R(unsigned tid,const dcmplx* in,dcmplx* out);

        void wait_for_workers(std::vector<std::thread>& workers);

        void hSTR_To_STSR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out);
        void hSTSR_To_STR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out);

};

*/




/*


// Attempt at using a third party thread pool - slow results :(

void worker_RT_To_RST(const double* in,dcmplx* out,\
        std::map<std::thread::id,FFTRVT*> fftblt_map);
void worker_STR_To_STSR(const dcmplx* in,dcmplx* out,\
        std::map<std::thread::id,HankelTransformR*> hankel_map);


class HankelFFTRRVTc 
{
    public:
        explicit HankelFFTRRVTc(const BesselRootGridR& gridR,const UniformGridT& gridT);
        ~HankelFFTRRVTc();
        HankelFFTRRVTc()=delete;
        HankelFFTRRVTc(const HankelFFTRRVTc& sp)=delete;
        HankelFFTRRVTc& operator=(const HankelFFTRRVTc& sp)=delete;

        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 
        bool setThreadPool(thread_pool* pool);

    private:
        int m_nr;
        int m_nt;
        int m_nst;
        int m_nThreads;

        std::vector<double> m_rr;
        std::vector<dcmplx> m_rs;
        std::vector<double> m_sr;
        std::vector<dcmplx> m_ss;

        BesselRootGridR* m_gridR;
        UniformGridT* m_gridT;
        FFTRVT* m_fftblt;
        HankelTransformR* m_hankel;
        std::map<std::thread::id,FFTRVT*> m_transformT;
        std::map<std::thread::id,HankelTransformR*> m_transformR;
        thread_pool* m_pool;

        void worker_RT_To_RST(unsigned tid);
        void worker_STR_To_STSR(unsigned tid);

};

*/



}



