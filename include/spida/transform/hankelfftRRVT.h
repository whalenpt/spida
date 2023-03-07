// Class for 2D functions of radius (R) and band-limited time (RVT)
#pragma once

#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include "spida/helper/constants.h"
#include "spida/transform/fftRVT.h"
#include "spida/transform/hankelR.h"

namespace spida{

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
        unsigned m_nr;
        unsigned m_nt;
        unsigned m_nst;
        unsigned m_nThreads;

        std::vector<double> m_rr; std::vector<double> m_sr;
        std::vector<dcmplx> m_rs; std::vector<dcmplx> m_ss;

        std::vector<std::thread> m_thread;
        std::vector<std::unique_ptr<FFTRVT>> m_transformT;
        std::vector<std::unique_ptr<HankelTransformR>> m_transformR;

        State m_state{State::Wait};
        unsigned m_th_count{0};
        std::vector<bool> m_ready;
        std::mutex m_mut;
        std::condition_variable m_cv;
        void setState(State state);

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

}