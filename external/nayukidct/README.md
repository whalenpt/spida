# nayukidct #

c++ fast cosine transform (FCT) library

# Codebase #

Library built primarily from code provided by Nayuki 
([details here](https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms)) with some added
tests, demos, and CMake integration.

# Dependencies #

C++11 or later, CMake

# Usage #

Fast discrete cosine transform for real inputs of even symmetry.
Useful in spectral methods utilizing Chebyshev interpolation grids,
and in image compression algorithms, among other uses. An FFT can often
be utilized instead of a FCT, but the typical complex-to-complex FFT
algorithm will have around a 4-fold increase in redundancy and extra
computation cost (in addition to being a less natural fit for the problem at hand).

Algorithms are implemented using a method proposed by Byeong Gi Lee
(A New Algorithm to Compute the Discrete Cosine Transform). This method 
implementes an FCT of type-II (Chebyshev root basis) for the forward transform 
and an FCT of type-III for the inverse transform 
([more info](https://en.wikipedia.org/wiki/Discrete_cosine_transform)) . The fast indicates
an *O(NlogN)* computation vs the "not fast" DCT *O(N^2)* computation, where *N* is the number of
points.

An example using the FCT for taking a derivative (of exp(4*x)) is shown below

```cpp
#include <nayukidct/FastDctLee.hpp>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>

int main()
{
    int N = 64;
    double pi = 3.141592653589793238462643383279502884;

    // Set up Chebyshev Root Grid on [a,b]
    std::vector<double> x(N);
    double dc = pi/static_cast<double>(N);
    for(auto j = 0; j < N; j++) x[j] = cos((j+0.5)*dc); 

    // Setup exponential input
    std::vector<double> in(N);
    for(auto i = 0; i < x.size(); i++)
        in[i] = exp(4.0*x[i]);

    std::vector<double> alpha(N);
    std::copy(in.cbegin(),in.cend(),alpha.begin());
    FastDctLee::transform(alpha);
    for(auto& item : alpha)
        item /= static_cast<double>(N/2.0);

    // Take spectral derivative
    std::vector<double> beta(N);
    beta[N-1] = 0.0;
    beta[N-2] =  2.0*(N-1)*alpha[N-1]; 
    for(int k = N-2; k > 0; k--)
        beta[k-1] = beta[k+1] + 2.0*k*alpha[k];

    FastDctLee::inverseTransform(beta);

    double diff = 0.0;
    for(auto i = 0; i < in.size(); i++)
        diff += abs(4.0*in[i] - beta[i]);
    diff /= N;
    assert(diff < 1.0e-10);
    std::cout << "The derivative of exp(4*x) was correctly found to be 4exp(4*x)" << std::endl;
    return 0;
}
```

# Installation #
From the github source with cmake
```bash

git clone https://github.com/whalenpt/nayukidct.git
cd nayukidct
cmake -S . -B build 
cd build
make -j4
cmake --install .
```

# License #
This project is licensed under the MIT License - see the [LICENSE](./LICENSE.txt) file for details.

# Contact # 
Patrick Whalen - whalenpt@gmail.com


