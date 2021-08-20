
#ifndef SPIDA_TRANSFORM_HANKELR_H_
#define SPIDA_TRANSFORM_HANKELR_H_ 

#include <vector>
#include <iostream>
#include "spida/transform/transformR.h"
#include "spida/grid/besselR.h" 

namespace spida{

class HankelTransformR;
void printHankel(const HankelTransformR& transform,std::ostream& os = std::cout);

class HankelTransformR : public TransformR
{
    public:
        explicit HankelTransformR(const BesselRootGridR& grid);
        ~HankelTransformR() {};
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
        int getNr() const {return m_nr;}
    private:
        int m_nr;
        double m_alpha;
        std::vector<double> m_Ymk; 
        void initDHT(const BesselRootGridR& grid);
};

}




#endif // End of SPIDA_TRANSFORM_HANKELR_H_



