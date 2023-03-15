#pragma once

#include <vector>
#include <iostream>
#include "spida/helper/constants.h"

namespace spida{

class BesselRootGridR;
class HankelTransformR 
{
    public:
        explicit HankelTransformR(const BesselRootGridR& grid);
        ~HankelTransformR() = default;
        HankelTransformR()=delete;
        HankelTransformR(const HankelTransformR& sp)=delete;
        HankelTransformR& operator=(const HankelTransformR& sp)=delete;

        // Read data input to transform
        void R_To_SR(const double* in,double* out);
        void SR_To_R(const double* in,double* out);
        void R_To_SR(const std::vector<double>& in,std::vector<double>& out) 
            {return R_To_SR(in.data(),out.data());}
        void SR_To_R(const std::vector<double>& in,std::vector<double>& out)
            { return SR_To_R(in.data(),out.data());}
 
        // Complex data input to transform
        void R_To_SR(const dcmplx* in,dcmplx* out);
        void SR_To_R(const dcmplx* in,dcmplx* out);
        void R_To_SR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
            { return R_To_SR(in.data(),out.data());}
        void SR_To_R(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
            { return SR_To_R(in.data(),out.data());}

        const std::vector<double>& getYmk() const {return m_Ymk;}
        unsigned getNr() const {return m_nr;}
    private:
        unsigned m_nr;
        double m_alpha;
        std::vector<double> m_Ymk; 
        std::vector<dcmplx> m_YmkC; 
        void initDHT(const BesselRootGridR& grid);
};

class HankelTransformRb 
{
    public:
        explicit HankelTransformRb(const BesselRootGridR& grid,unsigned threads=1);
        ~HankelTransformRb() = default;
        HankelTransformRb()=delete;
        HankelTransformRb(const HankelTransformRb& sp)=delete;
        HankelTransformRb& operator=(const HankelTransformRb& sp)=delete;

        // Read data input to transform
        void R_To_SR(const double* in,double* out);
        void SR_To_R(const double* in,double* out);
        void R_To_SR(const std::vector<double>& in,std::vector<double>& out) 
            {return R_To_SR(in.data(),out.data());}
        void SR_To_R(const std::vector<double>& in,std::vector<double>& out)
            { return SR_To_R(in.data(),out.data());}
 
        // Complex data input to transform
        void R_To_SR(const dcmplx* in,dcmplx* out);
        void SR_To_R(const dcmplx* in,dcmplx* out);
        void R_To_SR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) 
            { return R_To_SR(in.data(),out.data());}
        void SR_To_R(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
            { return SR_To_R(in.data(),out.data());}

        const std::vector<double>& getYmk() const {return m_Ymk;}
        unsigned getNr() const {return m_nr;}
    private:
        unsigned m_threads;
        unsigned m_nr;
        double m_alpha;
        std::vector<double> m_Ymk; 
        void initDHT(const BesselRootGridR& grid);
};

}