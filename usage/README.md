
# Usage #

First install the SPIDA library (see Installation instructions in the top-level README file).
Then, using CMake, in your CMakeLists.txt file all that is needed two lines of commands

1. find_package(SPIDA) # such that your program knows where SPIDA is installed 
2. target_link_libraries(target_name PRIVATE SPIDA::spida) 

An example of SPIDA usage is in the main.cpp file. To build this example using
cmake (after installing SPIDA), run the following from the command line
```bash
cmake -S . -B build
cd build
make
```
An executable called usage is created in the build directory. 
To run this executable:
```bash
./usage
```





