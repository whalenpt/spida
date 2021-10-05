
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
    if(sign_reverse){
        for(int i = in.size()-1; i >=0; i--)
            out[(in.size()-1)-i] = -in[i];
    } else{
        for(int i = in.size()-1; i >=0; i--)
            out[(in.size()-1)-i] = in[i];
    }
    for(auto i = 0; i < in.size(); i++)
        out[i+in.size()] = in[i];
}


void GridR::mirrorGrid(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,bool sign_reverse) const
{
    if(out.size() != 2*in.size())
        out.resize(2*in.size());
    if(sign_reverse){
        for(int i = in.size()-1; i >=0; i--)
            out[(in.size()-1)-i] = -in[i];
    } else{
        for(int i = in.size()-1; i >=0; i--)
            out[(in.size()-1)-i] = in[i];
    }
    for(auto i = 0; i < in.size(); i++)
        out[i+in.size()] = in[i];
}


}










