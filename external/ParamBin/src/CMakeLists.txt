
cmake_minimum_required(VERSION 3.19.1)

set(SRC_FILES
	parambin.cpp
	parambinstrings.cpp
)

set(INCL_DIR "${ParamBin_SOURCE_DIR}/include/ParamBin")
set(INCLUDE_FILES
    "${INCL_DIR}/parambin.hpp"
    "${INCL_DIR}/parambinstrings.h"
)

add_library(parambin ${SRC_FILES} ${INCLUDE_FILES})
add_library(ParamBin::parambin ALIAS parambin)
target_include_directories(parambin PRIVATE
	$<BUILD_INTERFACE:${ParamBin_SOURCE_DIR}/include/ParamBin>
	$<INSTALL_INTERFACE:include/ParamBin>
)
target_include_directories(parambin PUBLIC 
	$<BUILD_INTERFACE:${ParamBin_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

set_target_properties(parambin PROPERTIES
	VERSION ${ParamBin_VERSION}
	SOVERSION 0
)

install(TARGETS parambin
	EXPORT ParamBinTargets
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
)


install(DIRECTORY ${INCL_DIR}
	DESTINATION include
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	ParamBinConfigVersion.cmake
	VERSION ${ParamBin_VERSION}
	COMPATIBILITY AnyNewerVersion
)

install(EXPORT ParamBinTargets
	FILE ParamBinTargets.cmake
	NAMESPACE ParamBin::
	DESTINATION lib/cmake/ParamBin
)

configure_file(${ParamBin_SOURCE_DIR}/cmake/ParamBinConfig.cmake.in 
    ParamBinConfig.cmake @ONLY
)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/ParamBinConfig.cmake
	      ${CMAKE_CURRENT_BINARY_DIR}/ParamBinConfigVersion.cmake
	DESTINATION lib/cmake/ParamBin
)



