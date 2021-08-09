
#include <fftw3.h>
#include <cmath>
#include "spida/transform/fftwT.h"
#include "spida/grid/gridT.h" 
#include "spida/constants.h"

namespace spida{

PeriodicTransformT::PeriodicTransformT(UniformGridT& grid) 
{
    nt = grid.getNT();
    nst = grid.getNST();
    minI = grid.getMinI();
    maxI = grid.getMaxI();
    uFFT = std::unique_ptr<double[]>(new double[nt+2]);
    uFFTc = std::unique_ptr<dcmplx[]>(new dcmplx[nt]);
    for(int i = 0; i < nt+2; i++) uFFT[i] = 0.0;
    for(int i = 0; i < nt; i++) uFFTc[i] = 0.0;
    re_fftw_plan = fftw_plan_dft_r2c_1d(nt,uFFT.get(),reinterpret_cast<fftw_complex*>(uFFT.get()), FFTW_ESTIMATE);
    sp_fftw_plan = fftw_plan_dft_c2r_1d(nt,reinterpret_cast<fftw_complex*>(uFFT.get()), uFFT.get(), FFTW_ESTIMATE);
    re_fftw_plan_c = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(uFFTc.get()),reinterpret_cast<fftw_complex*>(uFFTc.get()),FFTW_BACKWARD, FFTW_ESTIMATE);
    sp_fftw_plan_c = fftw_plan_dft_1d(nt,reinterpret_cast<fftw_complex*>(uFFTc.get()),reinterpret_cast<fftw_complex*>(uFFTc.get()), FFTW_FORWARD, FFTW_ESTIMATE);
}


PeriodicTransformT::~PeriodicTransformT()
{
    fftw_destroy_plan(re_fftw_plan);
    fftw_destroy_plan(sp_fftw_plan);
    fftw_destroy_plan(re_fftw_plan_c);
    fftw_destroy_plan(sp_fftw_plan_c);

}

void PeriodicTransformT::T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out)
{
    // FFTW FORWARD - > +iwt transform def -> reverse time t -> -t
    for (int j = 0; j < nt; j++) 
      uFFT[j] = in[nt-j-1]; 
    uFFT[nt] = 0.0;
    uFFT[nt+1] = 0.0;
    fftw_execute(re_fftw_plan);
    for (int j = minI; j <= maxI; j++)  
        out[j-minI] = dcmplx(uFFT[2*j],uFFT[2*j+1]);
}

void PeriodicTransformT::ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out)
{
    // band limited minimum frequency 
    for(int j = 0; j < 2*minI; j++) 
        uFFT[j] = 0.0; 
    for(int j = minI; j <= maxI; j++){
        uFFT[2*j] = in[j-minI].real()/nt; 
        uFFT[2*j+1] = in[j-minI].imag()/nt; 
    }
    // band limited maximum frequency
    for(int j = 2*maxI+2; j < nt+2; j++) 
        uFFT[j] = 0.0; 
    fftw_execute(sp_fftw_plan);
    // FFTW BACKWARD - > -iwt transform def -> reverse time t -> -t
    for (int j = 0; j < nt; j++)  
        out[j] = uFFT[nt-j-1];
}

void PeriodicTransformT::T_To_ST_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    for (int j = 0; j < nt; j++) 
        uFFTc[j] = in[j];
    fftw_execute(re_fftw_plan_c);
    for (int j = minI; j <= maxI; j++)  
        out[j-minI] = uFFTc[j];
}

void PeriodicTransformT::ST_To_T_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out)
{
    for (int j = 0; j < minI; j++) 
        uFFTc[j] = dcmplx(0.0,0.0);
    for(int j = minI; j <= maxI; j++) 
        uFFTc[j] = in[j-minI]/static_cast<double>(nt); 
    for (int j = maxI+1; j < nt; j++)  
        uFFTc[j] = dcmplx(0.0,0.0);
    fftw_execute(sp_fftw_plan_c);
    for (int j = 0; j < nt; j++)  
        out[j] = uFFTc[j];
    
}



}


