
cmake_minimum_required(VERSION 3.19)

add_executable(interp_example interp_example.cpp)
target_link_libraries(interp_example PRIVATE 
    SPIDA::spida
    pw::pwutils
    ParamBin::parambin
)


