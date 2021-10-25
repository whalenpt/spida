
cmake_minimum_required(VERSION 3.19)

add_executable(cheb_der cheb_der.cpp)
target_link_libraries(cheb_der PUBLIC 
    nayukidct
)


