cmake_minimum_required(VERSION 3.19.1)
set(CMAKE_CXX_STANDARD 11)

project(NayukiDCT VERSION 0.0.1 LANGUAGES CXX)

option(DCT_TEST "Build NayukiDCT tests" ON)
option(DCT_DEMOS "Build NayukiDCT demos" ON)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    include(GNUInstallDirs)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

add_subdirectory(src)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    if(DCT_TEST)
        enable_testing()
        set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)")
        add_subdirectory(external/googletest)
        add_subdirectory(test)
    endif()
    if(DCT_DEMOS)
        add_subdirectory(demos)
    endif()
endif()


