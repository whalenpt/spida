
# Usage #

First install the ParamBin library (see Installation instructions in the top-level README file).
Then, using CMake, in your CMakeLists.txt file all that is needed two lines of commands

1. find_package(ParamBin) # such that your program knows where ParamBin is installed 
2. target_link_library(target_name PRIVATE ParamBin::parambin) 

If using the shared ParamBin::parambin library on Windows, (1) the installed libparam.dll file
must be copied into the UseParamBin build folder or; (2) the bin folder containing libparam.dll
file must be added to the environment path (e.g. C:\Program Files (x86)\ParamBin\bin can be added
to the system path).

A simple example of ParamBin usage is in the main.cpp file. To build this example using
cmake, run the following from the command line
```bash
cmake -S . -B build
cd build
make
```
An executable called use_parambin is created in the build directory. 
To run this executable:
```bash
./use_parambin
```





