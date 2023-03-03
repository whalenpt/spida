#!/bin/bash
rm -rf build && mkdir build && cd build
cmake -S .. -DSPIDA_TEST=1 \
            -DSPIDA_DEMOS=1 \
            -DSPIDA_STATIC=1 \
            -DCMAKE_EXPORT_COMPILE_COMMANDS=1
cmake --build . --parallel 4
