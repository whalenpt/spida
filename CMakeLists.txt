cmake_minimum_required(VERSION 3.28)

project(SPIDA
    VERSION 0.1.0
    LANGUAGES CXX
)

# -----------------------------
# Options
# -----------------------------
option(SPIDA_TEST "Build spida tests" OFF)
option(SPIDA_DEMOS "Build spida demos" OFF)
option(WITH_OPENBLAS "Use OpenBLAS" OFF)
option(SPIDA_COVERAGE "Enable code coverage instrumentation" OFF)
option(SPIDA_UPDATE_SUBMODULES "Update git submodules on configure" ON)
option(SPIDA_SETUP_GIT_HOOKS "Wire git core.hooksPath to .githooks on configure" ON)

set(KISSFFT_TOOLS OFF CACHE BOOL "" FORCE)
set(KISSFFT_PKGCONFIG OFF CACHE BOOL "" FORCE)
set(KISSFFT_TEST OFF CACHE BOOL "" FORCE)
set(KISSFFT_DATATYPE "double" CACHE STRING "" FORCE)

# Emit compile_commands.json for clangd / tooling.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
# Static build
# -----------------------------
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(KISSFFT_STATIC ON CACHE BOOL "" FORCE)

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
# Git hooks (dev convenience; only when developing spida itself)
# -----------------------------
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git"
        AND PROJECT_IS_TOP_LEVEL AND SPIDA_SETUP_GIT_HOOKS)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} config --local --get core.hooksPath
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE SPIDA_HOOKS_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT SPIDA_HOOKS_PATH STREQUAL ".githooks")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} config --local core.hooksPath .githooks
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            RESULT_VARIABLE SPIDA_HOOKS_RESULT
        )
        if(SPIDA_HOOKS_RESULT EQUAL 0)
            message(STATUS "Git hooks wired to .githooks/")
        else()
            message(WARNING "Failed to set core.hooksPath (${SPIDA_HOOKS_RESULT})")
        endif()
    endif()
endif()

# -----------------------------
# Dependencies
# -----------------------------
find_package(Threads REQUIRED)
find_package(Boost CONFIG QUIET)
find_package(spdlog CONFIG QUIET)

add_subdirectory(external/nayukidct)

find_package(kissfft CONFIG QUIET)
if(kissfft_FOUND)
    message(STATUS "Found kissfft package")
else()
    message(STATUS "Building kissfft from submodule")
    add_subdirectory(external/kissfft)
endif()

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
# clangd: copy compile_commands.json to source root (auto-detected there).
# Build-time copy so it works for any build dir (build/Release, build/RelWithDebInfo).
# -----------------------------
if(PROJECT_IS_TOP_LEVEL AND CMAKE_EXPORT_COMPILE_COMMANDS)
    add_custom_target(copy-compile-commands ALL
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_BINARY_DIR}/compile_commands.json"
                "${PROJECT_SOURCE_DIR}/compile_commands.json"
        DEPENDS "${CMAKE_BINARY_DIR}/compile_commands.json"
        COMMENT "Copying compile_commands.json to source root for clangd"
        VERBATIM)
endif()

# -----------------------------
# Tests / demos
# -----------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    if(SPIDA_TEST)
        enable_testing()
        find_package(GTest CONFIG QUIET)
        if(GTest_FOUND)
            message(STATUS "Found GTest package")
        else()
            message(STATUS "Building googletest from submodule")
            set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
            add_subdirectory(external/googletest)
        endif()
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
message(STATUS "Tests: ${SPIDA_TEST}")
message(STATUS "Demos: ${SPIDA_DEMOS}")
message(STATUS "Coverage: ${SPIDA_COVERAGE}")
