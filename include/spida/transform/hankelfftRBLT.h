
#ifndef SPIDA_TRANSFORM_HANKELFFTRBLT_H_
#define SPIDA_TRANSFORM_HANKELFFTRBLT_H_ 

#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <future>
#include <map>
#include "spida/helper/constants.h"

class thread_pool;

namespace spida{

class FFTBLT;
class HankelTransformR;
class HankelTransformRb;
class BesselRootGridR;
class UniformGridT;
void worker_RT_To_RST(const double* in,dcmplx* out,\
        std::map<std::thread::id,FFTBLT*> fftblt_map);
void worker_STR_To_STSR(const dcmplx* in,dcmplx* out,\
        std::map<std::thread::id,HankelTransformR*> hankel_map);

class HankelFFTRBLTb 
{
    public:
        explicit HankelFFTRBLTb(const BesselRootGridR& gridR,const UniformGridT& gridT);
        ~HankelFFTRBLTb();
        HankelFFTRBLTb()=delete;
        HankelFFTRBLTb(const HankelFFTRBLTb& sp)=delete;
        HankelFFTRBLTb& operator=(const HankelFFTRBLTb& sp)=delete;

        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 
        bool setThreadPool(thread_pool* pool);

//        void RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out); 
//        void RT_To_SRT(const std::vector<double>& in,std::vector<double>& out); 
//        void SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
//        void SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out);

    private:
        int m_nr;
        int m_nt;
        int m_nst;
        int m_nThreads;

        std::vector<double> m_uRT;
        std::vector<dcmplx> m_uRST;
        std::vector<double> m_uSRT;
        std::vector<dcmplx> m_uSRST;

        BesselRootGridR* m_gridR;
        UniformGridT* m_gridT;
        FFTBLT* m_fftblt;
        HankelTransformR* m_hankel;
        std::map<std::thread::id,FFTBLT*> m_transformT;
        std::map<std::thread::id,HankelTransformR*> m_transformR;
        thread_pool* m_pool;

        void worker_RT_To_RST(unsigned int tid);
        void worker_STR_To_STSR(unsigned int tid);


//        Semaphore m_task_limiter;

//        void worker_RT_To_RST(int tid);
//        void worker_RST_To_STR(int tid);
//        void worker_STR_To_STSR(int tid);
//        void worker_STSR_To_SRST(int tid);
//
//        void worker_SRST_To_STSR(int tid);
//        void worker_STSR_To_STR(int tid);
//        void worker_STR_To_RST(int tid);
//        void worker_RST_To_RT(int tid);
//
//        void worker_RT_To_SRT(int tid);
//        void worker_SRST_To_SRT(int tid);
};




class Semaphore{
    public:
        explicit Semaphore(unsigned int count) : m_count(count) {}
        ~Semaphore(){
            {
                std::unique_lock<std::mutex> lock(m_mut);
                m_count = 0;
            }
            m_cv.notify_all();
        }
        void lock() {
            std::unique_lock<std::mutex> lock(m_mut);
            m_cv.wait(lock, [this] {return m_count != 0;});
            m_count--;
        }

        void unlock() {
            std::unique_lock<std::mutex> lock(m_mut);
            m_count++;
            m_cv.notify_one();
        }
    private:
        std::mutex m_mut;
        std::condition_variable m_cv;
        unsigned int m_count;
};

template<typename T>
bool isReady(const std::future<T>& f) {
    if (f.valid()) { // otherwise you might get an exception (std::future_error: No associated state)
        return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    } else {
        return false;
    }
}

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

        std::vector<double> m_uRT;
        std::vector<dcmplx> m_uRST;
        std::vector<double> m_uSRT;
        std::vector<dcmplx> m_uSRST;

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


}


#endif // SPIDA_TRANSFORM_HANKELPERIODICRT_H_


