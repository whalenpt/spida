cmake_minimum_required(VERSION 3.21)

project(spida_test_package CXX)

# Matches the case of the config files spida actually installs
# (spidaConfig.cmake / spidaConfigVersion.cmake) — find_package(SPIDA) would
# not locate them on a case-sensitive filesystem.
find_package(spida REQUIRED)

add_executable(spida_test_package src/example.cpp)
target_link_libraries(spida_test_package PRIVATE SPIDA::spida)
