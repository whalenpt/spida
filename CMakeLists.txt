cmake_minimum_required(VERSION 3.19.1)
set(CMAKE_CXX_STANDARD 17)

project(Spida VERSION 0.0.1 LANGUAGES CXX)

option(SPIDA_STATIC "Build spida as static (ON) or shared library (OFF)" OFF)
option(SPIDA_TEST "Build spida tests" OFF)
option(SPIDA_EXAMPLES "Build spida examples" OFF)

set(PARAMBIN_TEST OFF CACHE BOOL "Build ParamBin tests")
set(PARAMBIN_EXAMPLES OFF CACHE BOOL "Build ParamBin examples")
set(PWUTILS_TEST OFF CACHE BOOL "Build pwutils tests")
set(PWUTILS_EXAMPLES OFF CACHE BOOL "Build pwutils examples")

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
	include(GNUInstallDirs)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
    endif()
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

if (NOT SPIDA_STATIC)
    set(BUILD_SHARED_LIBS ON)
    message(STATUS "Building shared spida library")
else()
    set(BUILD_SHARED_LIBS OFF)
    message(STATUS "Building static spida library")
    set(PWUTILS_STATIC ON CACHE BOOL "Build pwutils as static (ON) or shared library (OFF)")
    set(KISSFFT_STATIC ON CACHE BOOL "Build kissfft as static (ON) or shared library (OFF)")
    set(PARAMBIN_STATIC ON CACHE BOOL "Build ParamBin as static (ON) or shared library (OFF)")
endif()

if(NOT TARGET parambin)
    add_subdirectory(external/ParamBin)
endif()
if(NOT TARGET pwutils)
    add_subdirectory(external/pwutils)
endif()
if(NOT TARGET nayukidct)
    add_subdirectory(external/nayukidct)
endif()

set(JSON_DIR "external/json")
set(JSON_BuildTests OFF CACHE INTERNAL "")
# If you only include this third party in PRIVATE source files, you do not
# need to install it when your main project gets installed.
# set(JSON_Install OFF CACHE INTERNAL "")
add_subdirectory(${JSON_DIR})

#set(KISSFFT_PKGCONFIG OFF CACHE BOOL "Build pkg-config files")
set(KISSFFT_TOOLS OFF CACHE BOOL "Build kissfft command-line tools")
set(KISSFFT_PKGCONFIG OFF CACHE BOOL "Build pkg-config files")
set(KISSFFT_TEST OFF CACHE BOOL "Build kissfft tests")
set(KISSFFT_DIR "external/kissfft")
set(KISSFFT_DATATYPE "double" CACHE STRING "Principal datatype of kissfft: double, float (default), int16_t, int32_t, simd")
add_subdirectory(${KISSFFT_DIR})

add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if(SPIDA_TEST)
        enable_testing()
        set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)")
        add_subdirectory(external/googletest)
        add_subdirectory(test)
    endif()
    if(SPIDA_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()







