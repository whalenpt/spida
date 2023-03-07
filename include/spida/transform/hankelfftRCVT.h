#pragma once

#include <memory>
#include <thread>
#include <vector>
#include "spida/helper/constants.h"
#include "spida/transform/fftCVT.h"
#include "spida/transform/hankelR.h"

namespace spida{

class BesselRootGridR;
class UniformGridCVT;

class HankelFFTRCVT
{
    public:
        explicit HankelFFTRCVT(const BesselRootGridR& gridR,const UniformGridCVT& gridT,unsigned threads = 1);
        ~HankelFFTRCVT() = default;
        HankelFFTRCVT()=delete;
        HankelFFTRCVT(const HankelFFTRCVT& sp)=delete;
        HankelFFTRCVT& operator=(const HankelFFTRCVT& sp)=delete;

        void RT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 

        void RT_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void RST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void RST_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_RST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);

        void RT_To_SRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void SRT_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void SRT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
        void SRST_To_SRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);


    private:
        unsigned m_nr;
        unsigned m_nt;
        unsigned m_nst;
        unsigned m_threads;

        std::vector<dcmplx> m_rs;
        std::vector<dcmplx> m_ss;
        std::vector<dcmplx> m_rr;
        std::vector<dcmplx> m_sr;

        std::vector<std::unique_ptr<FFTCVT>> m_transformT;
        std::vector<std::unique_ptr<HankelTransformR>> m_transformR;

        void worker_T_To_ST(unsigned tid,const dcmplx* in,dcmplx* out);
        void worker_ST_To_T(unsigned tid,const dcmplx* in,dcmplx* out);
        void worker_R_To_SR(unsigned tid,const dcmplx* in,dcmplx* out);
        void worker_SR_To_R(unsigned tid,const dcmplx* in,dcmplx* out);
        void workerST_R_To_SR(unsigned tid,const dcmplx* in,dcmplx* out);
        void workerST_SR_To_R(unsigned tid,const dcmplx* in,dcmplx* out);
        void wait_for_workers(std::vector<std::thread>& workers);
        void hSTR_To_STSR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out);
        void hSTSR_To_STR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out);
};

}