// funcs.hpp
#pragma once

#include <cassert>
#include <cstddef>

namespace spida {

template <typename T> void transpose(const T* in, T* out, std::size_t nd1, std::size_t nd2)
{
    assert(out != in);
    for (std::size_t i = 0; i < nd1; i++)
        for (std::size_t j = 0; j < nd2; j++)
            out[j * nd1 + i] = in[i * nd2 + j];
}

} // namespace spida
