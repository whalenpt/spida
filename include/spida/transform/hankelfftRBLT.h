
#ifndef SPIDA_TRANSFORM_HANKELFFTRBLT_H_
#define SPIDA_TRANSFORM_HANKELFFTRBLT_H_ 

#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <future>
#include <map>
#include "spida/helper/constants.h"

namespace spida{

class FFTBLT;
class HankelTransformR;
class BesselRootGridR;
class UniformGridT;

template<typename T>
void transpose(const T* in,T* out,unsigned int nd1,unsigned int nd2)
{
    for(unsigned int i = 0; i < nd1; i++)
        for(unsigned int j = 0; j < nd2; j++)
            out[j*nd1+i] = in[i*nd2+j];
}

class HankelFFTRBLT 
{
    public:
        explicit HankelFFTRBLT(const BesselRootGridR& gridR,const UniformGridT& gridT,unsigned int threads = 1);
        ~HankelFFTRBLT();
        HankelFFTRBLT()=delete;
        HankelFFTRBLT(const HankelFFTRBLT& sp)=delete;
        HankelFFTRBLT& operator=(const HankelFFTRBLT& sp)=delete;

        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 

        void RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void RST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 

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

        std::vector<FFTBLT*> m_transformT;
        std::vector<HankelTransformR*> m_transformR;

        void worker_T_To_ST(unsigned int tid,const double* in,dcmplx* out);
        void worker_ST_To_T(unsigned int tid,const dcmplx* in,double* out);
        void worker_R_To_SR(unsigned int tid,const double* in,double* out);
        void worker_SR_To_R(unsigned int tid,const double* in,double* out);
        void workerCMP_R_To_SR(unsigned int tid,const dcmplx* in,dcmplx* out);
        void workerCMP_SR_To_R(unsigned int tid,const dcmplx* in,dcmplx* out);

        void wait_for_workers(std::vector<std::thread>& workers);

        void hSTR_To_STSR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out);
        void hSTSR_To_STR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out);

};


/*

class HankelFFTRBLT 
{
    public:
        explicit HankelFFTRBLT(const BesselRootGridR& gridR,const UniformGridT& gridT,int threads=1);
        ~HankelFFTRBLT();
        HankelFFTRBLT()=delete;
        HankelFFTRBLT(const HankelFFTRBLT& sp)=delete;
        HankelFFTRBLT& operator=(const HankelFFTRBLT& sp)=delete;
        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 

        void RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void RT_To_SRT(const std::vector<double>& in,std::vector<double>& out); 
        void SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out);

        enum class State{Wait,\
            RT_To_RST,\
            RST_To_STR,\
            STR_To_STSR,\
            STSR_To_SRST,\
            SRST_To_STSR,\
            STSR_To_STR,\
            STR_To_RST,\
            RST_To_RT,\
            RT_To_SRT,\
            SRST_To_SRT,\
            Done};

    private:
        int m_nr;
        int m_nt;
        int m_nst;
        int m_nThreads;

        std::vector<double> m_rr;
        std::vector<dcmplx> m_rs;
        std::vector<double> m_sr;
        std::vector<dcmplx> m_ss;

        std::vector<std::thread> m_thread;
        std::vector<FFTBLT*> m_transformT;
        std::vector<HankelTransformR*> m_transformR;

        State m_STATE;
        int m_THCOUNT;
        std::vector<bool> m_ready;
        bool m_processed;
        std::mutex m_mut;
        std::condition_variable m_cv;

        void ReadySTATE(State state);
        void ProcessedSTATE(State state);

        void worker_thread(int id);
        void worker_wait(int id);

        void worker_RT_To_RST(int tid);
        void worker_RST_To_STR(int tid);
        void worker_STR_To_STSR(int tid);
        void worker_STSR_To_SRST(int tid);

        void worker_SRST_To_STSR(int tid);
        void worker_STSR_To_STR(int tid);
        void worker_STR_To_RST(int tid);
        void worker_RST_To_RT(int tid);

        void worker_RT_To_SRT(int tid);
        void worker_SRST_To_SRT(int tid);
};
*/


/*


// Attempt at using a third party thread pool - slow results :(

void worker_RT_To_RST(const double* in,dcmplx* out,\
        std::map<std::thread::id,FFTBLT*> fftblt_map);
void worker_STR_To_STSR(const dcmplx* in,dcmplx* out,\
        std::map<std::thread::id,HankelTransformR*> hankel_map);


class HankelFFTRBLTc 
{
    public:
        explicit HankelFFTRBLTc(const BesselRootGridR& gridR,const UniformGridT& gridT);
        ~HankelFFTRBLTc();
        HankelFFTRBLTc()=delete;
        HankelFFTRBLTc(const HankelFFTRBLTc& sp)=delete;
        HankelFFTRBLTc& operator=(const HankelFFTRBLTc& sp)=delete;

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
        FFTBLT* m_fftblt;
        HankelTransformR* m_hankel;
        std::map<std::thread::id,FFTBLT*> m_transformT;
        std::map<std::thread::id,HankelTransformR*> m_transformR;
        thread_pool* m_pool;

        void worker_RT_To_RST(unsigned int tid);
        void worker_STR_To_STSR(unsigned int tid);

};

*/



}


#endif // SPIDA_TRANSFORM_HANKELPERIODICRT_H_


