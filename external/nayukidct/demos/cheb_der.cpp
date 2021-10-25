/*
MIT License

Copyright (c) 2021 Patrick Whalen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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



