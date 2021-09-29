

cmake_minimum_required(VERSION 3.19.1)

set(SRC_FILES
    SpidaBLT.cpp
    SpidaRBLT.cpp
    SpidaRVX.cpp
    SpidaCVX.cpp
    SpidaCHEBX.cpp
    ChebInterpX.cpp
	helper/interp.cpp
    grid/chebX.cpp
    grid/uniformX.cpp
    grid/uniformRVX.cpp
    grid/uniformCVX.cpp
	grid/uniformT.cpp
	grid/besselR.cpp
    grid/grid.cpp
	report/reportcenter.cpp
	report/reporthandler.cpp
	shape/shapeT.cpp
	shape/shapeR.cpp
	rkstiff/solver.cpp
	rkstiff/ETDAS.cpp
	rkstiff/ETDCS.cpp
    transform/fftBLT.cpp
    transform/fftCVX.cpp
    transform/fftRVX.cpp
	transform/chebX.cpp
	transform/hankelR.cpp
    transform/hankelfftRBLT.cpp
)

set(INCL_DIR "${Spida_SOURCE_DIR}/include/spida")
set(SPIDA_INCLUDES
   "${INCL_DIR}/SpidaRVX.h"
   "${INCL_DIR}/SpidaCVX.h"
   "${INCL_DIR}/SpidaCHEBX.h"
   "${INCL_DIR}/ChebInterpX.h"
   "${INCL_DIR}/SpidaBLT.h"
   "${INCL_DIR}/SpidaRBLT.h"
   "${INCL_DIR}/helper/constants.h"
   "${INCL_DIR}/helper/interp.h"
   "${INCL_DIR}/grid/grid.h"
   "${INCL_DIR}/grid/gridX.h"
   "${INCL_DIR}/grid/gridT.h"
   "${INCL_DIR}/grid/gridR.h"
   "${INCL_DIR}/grid/uniformX.h"
   "${INCL_DIR}/grid/uniformRVX.h"
   "${INCL_DIR}/grid/uniformCVX.h"
   "${INCL_DIR}/grid/uniformT.h"
   "${INCL_DIR}/grid/chebX.h"
   "${INCL_DIR}/grid/besselR.h"
   "${INCL_DIR}/propagator/propagator.h"
   "${INCL_DIR}/report/reportcenter.h"
   "${INCL_DIR}/report/reporthandler.h"
   "${INCL_DIR}/shape/shape.h"
   "${INCL_DIR}/shape/shapeT.h"
   "${INCL_DIR}/shape/shapeR.h"
   "${INCL_DIR}/rkstiff/solver.h"
   "${INCL_DIR}/rkstiff/ETDAS.h"
   "${INCL_DIR}/rkstiff/ETDCS.h"
   "${INCL_DIR}/transform/fftBLT.h"
   "${INCL_DIR}/transform/chebX.h"
   "${INCL_DIR}/transform/fftCVX.h"
   "${INCL_DIR}/transform/fftRVX.h"
   "${INCL_DIR}/transform/hankelR.h"
   "${INCL_DIR}/transform/hankelfftRBLT.h"
)

set(BOOST_DIR "${Spida_SOURCE_DIR}/external/boost_1_76_0")
set(BOOST_INCLUDES 
	"${BOOST_DIR}/boost/math/special_functions/airy.hpp"
	"${BOOST_DIR}/boost/math/special_functions/bessel.hpp"
)

set(INCLUDE_FILES
	${BOOST_INCLUDES}
	${SPIDA_INCLUDES}
)

add_library(spida ${SRC_FILES} ${INCLUDE_FILES})
add_library(SPIDA::spida ALIAS spida)

target_include_directories(spida PUBLIC 
    $<BUILD_INTERFACE:${Spida_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
	PRIVATE
	$<BUILD_INTERFACE:${BOOST_DIR}>
	$<INSTALL_INTERFACE:include/boost>
)

target_include_directories(spida PRIVATE
	$<BUILD_INTERFACE:${BOOST_DIR}>
	$<INSTALL_INTERFACE:include/boost>
    $<BUILD_INTERFACE:${Spida_SOURCE_DIR}/include/spida>
	$<INSTALL_INTERFACE:include/spida>
)

target_link_libraries(spida
    PRIVATE ParamBin::parambin
    PRIVATE pw::pwutils
    PUBLIC kissfft
    PUBLIC nayukidct
    #    PUBLIC gsl
)

set_target_properties(spida PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION 0
)

install(TARGETS spida
	EXPORT SpidaTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
)

install(DIRECTORY ${INCL_DIR}
    DESTINATION include
)

install(DIRECTORY "${BOOST_DIR}/boost"
	DESTINATION include
	PATTERN ".DS_Store" EXCLUDE
)

#install(DIRECTORY "${JSON_DIR}/include"
#	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#	PATTERN ".DS_Store" EXCLUDE
#)

install(EXPORT SpidaTargets
	FILE SpidaTargets.cmake
	NAMESPACE SPIDA::
	DESTINATION lib/cmake/spida
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	SpidaConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion
)
configure_file(${Spida_SOURCE_DIR}/cmake/SpidaConfig.cmake.in
    SpidaConfig.cmake @ONLY
)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfigVersion.cmake
    DESTINATION lib/cmake/spida
)







