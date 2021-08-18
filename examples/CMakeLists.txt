
cmake_minimum_required(VERSION 3.19)

add_executable(interp_ex interp_example.cpp)
add_executable(grid_ex grid_example.cpp)
add_executable(hankel_ex hankel_example.cpp)

target_link_libraries(interp_ex PUBLIC
    SPIDA::spida
    pw::pwutils
    ParamBin::parambin
)

target_link_libraries(grid_ex PUBLIC SPIDA::spida)
target_link_libraries(hankel_ex PUBLIC SPIDA::spida)



