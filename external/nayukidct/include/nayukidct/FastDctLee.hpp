/* 
 * Fast discrete cosine transform algorithms (C++)
 * 
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#pragma once

#include <cstddef>
#include <vector>

enum class dct_kind{DCTII,DCTIII};

class FastDctPlan
{
    public:
        FastDctPlan(std::vector<double>& in,std::vector<double>& out,dct_kind kind);
        void execute();
    private:
        std::size_t m_sz;
	    dct_kind m_type;
        std::vector<double>& m_in;
        std::vector<double>& m_out;
	    std::vector<double> m_temp;
	    bool m_inplace;
};


namespace FastDctLee {
    constexpr auto PI = 3.141592653589793238462643383279502884197;
	
	void transform(std::vector<double> &vec);
	
	void transform(double vec[], std::size_t len);
	
	void inverseTransform(std::vector<double> &vec);
	
	void inverseTransform(double vec[], std::size_t len);
	
}


