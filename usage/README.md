
# Usage #

First install the SPIDA library (see Installation instructions in the top-level README file).
Then, create a CMakeLists.txt file that includes the two lines of commands

1. find_package(SPIDA) # such that your program knows where SPIDA is installed 
2. target_link_libraries(target_name PRIVATE SPIDA::spida) 

An example of SPIDA usage is provided here with the given CMakeLists.txt file and
main.cpp source code. To build this example, make sure your in the usage directory
and have already installed SPIDA, then use the following terminal command
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





