cmake_minimum_required(VERSION 3.19)

set(SRC_FILES
    NaiveDct.cpp
    FastDct8.cpp
    FastDctLee.cpp
    FastDctFft.cpp
    FftRealPair.cpp
)

set(INCL_DIR "${NayukiDCT_SOURCE_DIR}/include/nayukidct")
set(INCLUDE_FILES
    "${INCL_DIR}/NaiveDct.hpp"
    "${INCL_DIR}/FastDct8.hpp"
    "${INCL_DIR}/FastDctLee.hpp"
    "${INCL_DIR}/FastDctFft.hpp"
    "${INCL_DIR}/FftRealPair.hpp"
)

add_library(nayukidct ${SRC_FILES} ${INCLUDE_FILES})
target_include_directories(nayukidct PRIVATE
    $<BUILD_INTERFACE:${NayukiDCT_SOURCE_DIR}/include/nayukidct>
    $<INSTALL_INTERFACE:include/nayukidct>
)

target_include_directories(nayukidct PUBLIC
    $<BUILD_INTERFACE:${NayukiDCT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(nayukidct PROPERTIES
	VERSION ${NayukiDCT_VERSION}
	SOVERSION 0
)

install(TARGETS nayukidct
    EXPORT NayukiDCTTargets
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	RUNTIME DESTINATION bin
)

install(DIRECTORY ${INCL_DIR}
	DESTINATION include
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
	NayukiDCTConfigVersion.cmake
	VERSION ${NayukiDCT_VERSION}
	COMPATIBILITY AnyNewerVersion
)

install(EXPORT NayukiDCTTargets
	FILE NayukiDCTTargets.cmake
	NAMESPACE pw::
	DESTINATION lib/cmake/nayukidct
)
configure_file(${NayukiDCT_SOURCE_DIR}/cmake/NayukiDCTConfig.cmake.in
    NayukiDCTConfig.cmake @ONLY
)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/NayukiDCTConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/NayukiDCTConfigVersion.cmake
    DESTINATION lib/cmake/nayukidct
)







