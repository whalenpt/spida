# ParamBin #

c++ class that acts as a tree structure for containing 
named-parameters of various data types.

# Dependencies #

ParamBin only requires standard c++ libraries (STL). It is 
installed with c++17 as a default and may need adjustments
to work with earlier c++ versions.

# Usage #

ParamBin class is designed for easy storage and access to parameters
that fit into a tree structure. 

```cpp
#include <ParamBin/parambin.hpp>
#include <string>
int main()
{
	ParamBin bin;
	bin.set("var1",1.3); // set a float
	bin.set("var2",10); // set an int
	bin.set("var3","string"); // set a string
	bin.setBool("var4",true); // set a bool

	ParamBin subbin;
	subbin.set("A",5.3e8);
	std::vector<std::string> str_vals;
	str_vals.push_back("first val");
	str_vals.push_back("second val");
	subbin.set("strvec",str_vals); // set a vector (of strings)
	bin.set("SUBBIN",subbin); // set a subbin of bin

	std::cout << bin << std::endl; // print the bin

	double var1 = bin.getDbl("var1"); // access a double
	int var2 = bin.getInt("var2"); // access an int
	std::string var3 = bin.getStr("var3"); // access a string
	ParamBin subbin_access = bin.getBin("SUBBIN"); // access a subbin
	std::vector<std::string> strvec_access = bin.getBin("SUBBIN").getStrVec("strvec"); // access a vector
}
```

# Installation #
From the github source with cmake
```bash
git clone https://github.com/whalenpt/ParamBin.git
cd ParamBin
cmake -S . -B build (-DPARAMBIN_STATIC=1 optional argument -> build static library)
cd build
make -j4
cmake --install .
```
A program using ParamBin will need to link to the parambin library
(or its static counterpart). See the usage folder
[README](./usage/README.md) for an example.

# License #
This project is licensed under the MIT License - see the [LICENSE](./LICENSE.txt) file for details

# Contact # 
Patrick Whalen - whalenpt@gmail.com

