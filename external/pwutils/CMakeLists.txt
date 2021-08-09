cmake_minimum_required(VERSION 3.19.1)
set(CMAKE_CXX_STANDARD 17)

project(pwutils VERSION 0.0.1 LANGUAGES CXX)

option(PWUTILS_STATIC "Build pwutils as static (ON) or shared library (OFF)" OFF)
option(PWUTILS_TEST "Build pwutils tests" OFF)
option(PWUTILS_EXAMPLES "Build pwutils examples" OFF)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
	include(GNUInstallDirs)
endif()

if (NOT PWUTILS_STATIC)
    set(BUILD_SHARED_LIBS ON)
    message(STATUS "Building shared pwutils library")
else()
    set(BUILD_SHARED_LIBS OFF)
    message(STATUS "Building static pwutils library")
endif()

add_subdirectory(src)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    if(PWUTILS_TEST)
        enable_testing()
        set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)")
        add_subdirectory(external/googletest)
        add_subdirectory(test)
    endif()
    if(PWUTILS_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()








