cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_STANDARD 17)

project(use_parambin)
find_package(ParamBin)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ParamBin::parambin) 


