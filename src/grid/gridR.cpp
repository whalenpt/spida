
#include "spida/grid/gridR.h"
#include <vector>
#include <iostream>

namespace spida{

std::vector<double> GridR::mirrorGrid(const std::vector<double>& in,bool sign_reverse) const
{
    std::vector<double> out(2*in.size());
    mirrorGrid(in,out,sign_reverse);
    return out;
}

std::vector<dcmplx> GridR::mirrorGrid(const std::vector<dcmplx>& in,bool sign_reverse) const
{
    std::vector<dcmplx> out(2*in.size());
    mirrorGrid(in,out,sign_reverse);
    return out;
}

void GridR::mirrorGrid(const std::vector<double>& in,std::vector<double>& out,bool sign_reverse) const
{
    if(out.size() != 2*in.size())
        out.resize(2*in.size());
    mirrorGrid(in.data(),out.data(),sign_reverse);
}


void GridR::mirrorGrid(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,bool sign_reverse) const
{
    if(out.size() != 2*in.size())
        out.resize(2*in.size());
    mirrorGrid(in.data(),out.data(),sign_reverse);
}


void GridR::mirrorGrid(const double* in,double* out,bool sign_reverse) const
{
    if(sign_reverse){
        for(int i = m_nr-1; i >=0; i--)
            out[(m_nr-1)-i] = -in[i];
    } else{
        for(int i = m_nr-1; i >=0; i--)
            out[(m_nr-1)-i] = in[i];
    }
    for(unsigned i = 0; i < m_nr; i++)
        out[i+m_nr] = in[i];
}


void GridR::mirrorGrid(const dcmplx* in,dcmplx* out,bool sign_reverse) const
{
    if(sign_reverse){
        for(int i = m_nr-1; i >=0; i--)
            out[(m_nr-1)-i] = -in[i];
    } else{
        for(int i = m_nr-1; i >=0; i--)
            out[(m_nr-1)-i] = in[i];
    }
    for(unsigned i = 0; i < m_nr; i++)
        out[i+m_nr] = in[i];
}

}