cmake_minimum_required(VERSION 3.19.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(SPIDA VERSION 0.0.1 LANGUAGES CXX)

option(SPIDA_STATIC "Build spida as static (ON) or shared library (OFF)" ON)
option(SPIDA_TEST "Build spida tests" OFF)
option(SPIDA_DEMOS "Build spida demos" OFF)

set(PWUTILS_TEST OFF CACHE BOOL "Build pwutils tests")
set(PWUTILS_EXAMPLES OFF CACHE BOOL "Build pwutils examples")

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}")
        endif()
    endif()
endif()


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
    set(PWUTILS_STATIC OFF CACHE BOOL "Build pwutils as static (ON) or shared library (OFF)")
    set(KISSFFT_STATIC OFF CACHE BOOL "Build kissfft as static (ON) or shared library (OFF)")
    set(PWUTILS_STATIC OFF CACHE BOOL "Build pwutils as static (ON) or shared library (OFF)")
    message(STATUS "Building shared spida library")
else()
    set(BUILD_SHARED_LIBS OFF)
    message(STATUS "Building static spida library")
    set(PWUTILS_STATIC ON CACHE BOOL "Build pwutils as static (ON) or shared library (OFF)")
    set(KISSFFT_STATIC ON CACHE BOOL "Build kissfft as static (ON) or shared library (OFF)")
    set(PWUTILS_STATIC ON CACHE BOOL "Build pwutils as static (ON) or shared library (OFF)")
endif()

if(NOT TARGET pwutils)
    add_subdirectory(external/pwutils)
endif()
if(NOT TARGET nayukidct)
    add_subdirectory(external/nayukidct)
endif()

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
    if(SPIDA_DEMOS)
        add_subdirectory(demos)
    endif()
endif()







