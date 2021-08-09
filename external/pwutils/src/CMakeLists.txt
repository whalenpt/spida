
cmake_minimum_required(VERSION 3.19)

set(SRC_FILES
	pwstats.cpp
	pwthreads.cpp
	pwstrings.cpp
	pwmath.cpp
	report/dat.cpp
	report/json.cpp
	report/basereport.cpp
	report/basetrack.cpp
	report/reporthelper.cpp
	report/pwdir.cpp
	report/dataio.cpp
)

set(INCL_DIR "${pwutils_SOURCE_DIR}/include/pwutils")
set(INCLUDE_FILES
   "${INCL_DIR}/pwstats.h"
   "${INCL_DIR}/pwthreads.h"
   "${INCL_DIR}/pwexcept.h"
   "${INCL_DIR}/pwconstants.h"
   "${INCL_DIR}/pwstrings.h"
   "${INCL_DIR}/pwmath.hpp"
   "${INCL_DIR}/report/pwdir.h"
   "${INCL_DIR}/report/basereport.h"
   "${INCL_DIR}/report/basetrack.hpp"
   "${INCL_DIR}/report/basedata.hpp"
   "${INCL_DIR}/report/dat.hpp"
   "${INCL_DIR}/report/json.hpp"
   "${INCL_DIR}/report/reporthelper.h"
   "${INCL_DIR}/report/dataio.hpp"
)

add_library(pwutils ${SRC_FILES} ${INCLUDE_FILES})
add_library(pw::pwutils ALIAS pwutils)

target_include_directories(pwutils PRIVATE
	$<BUILD_INTERFACE:${pwutils_SOURCE_DIR}/include/pwutils>
	$<INSTALL_INTERFACE:include/pwutils>
)

target_include_directories(pwutils PUBLIC 
	$<BUILD_INTERFACE:${pwutils_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

set_target_properties(pwutils PROPERTIES
	VERSION ${pwutils_VERSION}
	SOVERSION 0
)

install(TARGETS pwutils
	EXPORT PWUtilsTargets
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	RUNTIME DESTINATION bin
)

install(DIRECTORY ${INCL_DIR}
	DESTINATION include
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
	PWUtilsConfigVersion.cmake
	VERSION ${pwutils_VERSION}
	COMPATIBILITY AnyNewerVersion
)

install(EXPORT PWUtilsTargets
	FILE PWUtilsTargets.cmake
	NAMESPACE pw::
	DESTINATION lib/cmake/pwutils
)
configure_file(${pwutils_SOURCE_DIR}/cmake/PWUtilsConfig.cmake.in
    PWUtilsConfig.cmake @ONLY
)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/PWUtilsConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/PWUtilsConfigVersion.cmake
    DESTINATION lib/cmake/pwutils
)



