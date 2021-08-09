# pwutils #

C++ helper functions and classes for repetitive tasks

# Dependencies #

The pwutils library only depends on standard c++ libraries (STL). 
Its built with c++17 by default and may need adjustments
to work with earlier c++ versions.

Testing is done with the googletest library using the
cmake compile flag -DPWUTILS_TEST::BOOL=ON

# Usage #

pwutils functions and classes are used for commonly
repeated tasks. For example upper casing or
lower casing a string:

```cpp
#include <pwutils/pwstrings.h>
#include <string>
#include <iostream>
int main()
{
    std::string str("UPPER and lower CaSes");
    // This will convert str to all lower case letters
	std::cout << pw::stringLowerCase(str) << std::endl;
    // This will convert str to all upper case letters
	std::cout << pw::stringUpperCase(str) << std::endl;
	return 0;
}
```

# Installation #
From the github source with cmake
```bash
git clone https://github.com/whalenpt/pwutils.git
cd pwutils
cmake -S . -B build
cd build
make -j4
cmake --install .
```
For a static build use
```bash
cmake -S . -B build -DPWUTILS_STATIC::BOOL=ON
```
Tests are built using the cmake compiler option
```bash
cmake -S . -B build -DPWUTILS_TEST::BOOL=ON
```
Examples are built using the cmake compiler option
```bash
cmake -S . -B build -DPWUTILS_EXAMPLES::BOOL=ON
```
A program using pwutils will need to link to the pwutils library 

# License #
This project is licensed under the MIT License - see the [LICENSE](./LICENSE.txt) file for details

# Contact # 
Patrick Whalen - whalenpt@gmail.com

