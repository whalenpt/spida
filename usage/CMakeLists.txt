cmake_minimum_required(VERSION 3.19)

project(usage)
# Lowercase to match what's actually installed (spidaConfig.cmake) —
# find_package(SPIDA) doesn't find it on a case-sensitive filesystem.
find_package(spida REQUIRED)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE SPIDA::spida)


