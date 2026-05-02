cmake_minimum_required(VERSION 3.20)

project(SPIDA
    VERSION 0.1.0
    LANGUAGES CXX
)

# -----------------------------
# Options
# -----------------------------
option(SPIDA_STATIC "Build spida as static library" ON)
option(SPIDA_TEST "Build spida tests" OFF)
option(SPIDA_DEMOS "Build spida demos" OFF)
option(WITH_OPENBLAS "Use OpenBLAS" OFF)
option(SPIDA_COVERAGE "Enable code coverage instrumentation" OFF)
option(SPIDA_UPDATE_SUBMODULES "Update git submodules on configure" ON)

set(PWUTILS_TEST OFF CACHE BOOL "Build pwutils tests")
set(PWUTILS_EXAMPLES OFF CACHE BOOL "Build pwutils examples")
set(KISSFFT_TOOLS OFF CACHE BOOL "" FORCE)
set(KISSFFT_PKGCONFIG OFF CACHE BOOL "" FORCE)
set(KISSFFT_TEST OFF CACHE BOOL "" FORCE)
set(KISSFFT_DATATYPE "double" CACHE STRING "" FORCE)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# -----------------------------
# Policies
# -----------------------------
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# -----------------------------
# Build type default
# -----------------------------
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# -----------------------------
# Git submodules
# -----------------------------
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git" AND SPIDA_UPDATE_SUBMODULES)
    message(STATUS "Updating git submodules")

    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_RESULT
    )

    if(NOT GIT_RESULT EQUAL 0)
        message(FATAL_ERROR "git submodule update failed: ${GIT_RESULT}")
    endif()
endif()

# -----------------------------
# Dependencies
# -----------------------------
find_package(Boost CONFIG QUIET)
find_package(pwutils CONFIG QUIET)

if(NOT pwutils_FOUND)
    message(STATUS "Building pwutils from submodule")
    add_subdirectory(external/pwutils)
else()
    message(STATUS "Found pwutils package")
endif()

add_subdirectory(external/nayukidct)

add_subdirectory(external/kissfft)

# -----------------------------
# OpenBLAS (optional)
# -----------------------------
if(WITH_OPENBLAS)
    find_package(OpenBLAS QUIET)

    if(OpenBLAS_FOUND)
        message(STATUS "OpenBLAS found")

        add_library(OpenBLAS::interface INTERFACE IMPORTED)
        target_compile_definitions(OpenBLAS::interface INTERFACE HAVE_OPENBLAS=1)
        target_include_directories(OpenBLAS::interface INTERFACE ${OpenBLAS_INCLUDE_DIRS})
        target_link_libraries(OpenBLAS::interface INTERFACE ${OpenBLAS_LIBRARIES})
    endif()
endif()

# -----------------------------
# Coverage (FLAGS ONLY — IMPORTANT FIX)
# -----------------------------
if(SPIDA_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(STATUS "Code coverage enabled")

    add_compile_options(--coverage -O0 -g -fprofile-update=atomic)
    add_link_options(--coverage)
endif()

# -----------------------------
# Main library
# -----------------------------
add_subdirectory(src)

# -----------------------------
# Tests / demos
# -----------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    if(SPIDA_TEST)
        enable_testing()
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
        add_subdirectory(external/googletest)
        add_subdirectory(test)
    endif()

    if(SPIDA_DEMOS)
        add_subdirectory(demos)
    endif()

endif()

# -----------------------------
# Summary
# -----------------------------
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Static build: ${SPIDA_STATIC}")
message(STATUS "Tests: ${SPIDA_TEST}")
message(STATUS "Demos: ${SPIDA_DEMOS}")
message(STATUS "Coverage: ${SPIDA_COVERAGE}")
