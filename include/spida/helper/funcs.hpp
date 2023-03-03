//funcs.hpp
#pragma once

namespace spida{

template<typename T>
void transpose(const T* in,T* out,unsigned nd1,unsigned nd2)
{
    for(unsigned i = 0; i < nd1; i++)
        for(unsigned j = 0; j < nd2; j++)
            out[j*nd1+i] = in[i*nd2+j];
}

}