
cmake_minimum_required(VERSION 3.19)
project(json11 VERSION 1.0.0 LANGUAGES CXX)

option(JSON11_BUILD_TESTS "Build unit tests" OFF)
option(JSON11_ENABLE_DR1467_CANARY "Enable canary test for DR 1467" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
	include(GNUInstallDirs)
    include(CheckCXXCompilerFlag)
    set(_possible_warnings_flags /W4 /WX -Wextra -Werror)
    foreach(_warning_flag ${_possible_warnings_flags})
      unset(_flag_supported)
      CHECK_CXX_COMPILER_FLAG(${_warning_flag} _flag_supported)
      if(${_flag_supported})
        target_compile_options(json11 PRIVATE ${_warning_flag})
      endif()
    endforeach()
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()


add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if (JSON11_BUILD_TESTS)
        enable_testing()
        # enable test for DR1467, described here: https://llvm.org/bugs/show_bug.cgi?id=23812
        if(JSON11_ENABLE_DR1467_CANARY)
          add_definitions(-D JSON11_ENABLE_DR1467_CANARY=1)
        else()
          add_definitions(-D JSON11_ENABLE_DR1467_CANARY=0)
        endif()
        add_subdirectory(test)
    endif()
endif()




