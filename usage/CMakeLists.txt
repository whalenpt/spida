cmake_minimum_required(VERSION 3.19)

project(usage)
find_package(SPIDA)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE SPIDA::spida)


