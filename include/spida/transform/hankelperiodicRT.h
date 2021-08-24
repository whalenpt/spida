
#ifndef SPIDA_TRANSFORM_HANKELPERIODICRT_H_
#define SPIDA_TRANSFORM_HANKELPERIODICRT_H_ 

#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>
#include "spida/constants.h"
#include "spida/transform/transformRT.h"

namespace spida{

class PeriodicTransformT;
class HankelTransformR;
class BesselRootGridR;
class UniformGridT;

class HankelPeriodicTransformRT : public TransformsRT
{
    public:
        explicit HankelPeriodicTransformRT(const BesselRootGridR& gridR,const UniformGridT& gridT,int threads=1);
        ~HankelPeriodicTransformRT();
        HankelPeriodicTransformRT()=delete;
        HankelPeriodicTransformRT(const HankelPeriodicTransformRT& sp)=delete;
        HankelPeriodicTransformRT& operator=(const HankelPeriodicTransformRT& sp)=delete;
        void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 

        void RT_To_RST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void RT_To_SRT(const std::vector<double>& in,std::vector<double>& out); 
        void SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<double>& out);
        enum class State{Wait,RT_To_RST,RST_To_SRST,SRST_To_RST,RST_To_RT,RT_To_SRT,SRST_To_SRT,Done};

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
        std::vector<PeriodicTransformT*> m_transformT;
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
        void worker_RST_To_SRST(int tid);
        void worker_SRST_To_RST(int tid);
        void worker_RST_To_RT(int tid);

        void worker_RT_To_SRT(int tid);
        void worker_SRST_To_SRT(int tid);

};

}


#endif // SPIDA_TRANSFORM_HANKELPERIODICRT_H_


