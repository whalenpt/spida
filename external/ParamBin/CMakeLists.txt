cmake_minimum_required(VERSION 3.19.1)

project(ParamBin VERSION 0.0.3 DESCRIPTION "Container for parameters")
set(CMAKE_CXX_STANDARD 17)

option(PARAMBIN_TEST "Build ParamBin tests" ON)
option(PARAMBIN_EXAMPLES "Build ParamBin examples" ON)
option(PARAMBIN_STATIC "Build ParamBin as static (ON) or shared library (OFF)" OFF)

if (NOT PARAMBIN_STATIC)
    set(BUILD_SHARED_LIBS ON)
    message(STATUS "Building shared ParamBin library")
else()
    set(BUILD_SHARED_LIBS OFF)
    message(STATUS "Building static ParamBin library")
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(GNUInstallDirs)
    set(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# source files
add_subdirectory(src)

# only build testing and examples if this is at the project-level 
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if(PARAMBIN_TEST)
        enable_testing()
        set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)")
        add_subdirectory(external/googletest)
        add_subdirectory(test)
    endif()
    # executable examples 
    if(PARAMBIN_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()







