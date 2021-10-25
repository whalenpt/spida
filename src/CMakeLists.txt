

cmake_minimum_required(VERSION 3.19.1)

set(SRC_FILES
    SpidaCVT.cpp
    SpidaRVT.cpp
    SpidaRCVT.cpp
    SpidaRRVT.cpp
    SpidaCVX.cpp
    SpidaRVX.cpp
    SpidaCHEBX.cpp
    SpidaR.cpp
    ChebInterpX.cpp
	helper/interp.cpp
    grid/grid.cpp
    grid/gridR.cpp
    grid/chebX.cpp
    grid/uniformX.cpp
    grid/uniformCVX.cpp
    grid/uniformRVX.cpp
	grid/uniformT.cpp
    grid/uniformCVT.cpp
    grid/uniformRVT.cpp
	grid/besselR.cpp
	propagator/propagator.cpp
	propagator/reporthandler.cpp
	shape/shapeT.cpp
	shape/shapeR.cpp
	rkstiff/solver.cpp
	rkstiff/ETDAS.cpp
	rkstiff/ETDCS.cpp
    rkstiff/IFAS.cpp
    transform/fftCVT.cpp
    transform/fftRVT.cpp
    transform/fftCVX.cpp
    transform/fftRVX.cpp
	transform/chebX.cpp
	transform/hankelR.cpp
    transform/hankelfftRCVT.cpp
    transform/hankelfftRRVT.cpp
)

set(INCL_DIR "${SPIDA_SOURCE_DIR}/include/spida")
set(SPIDA_INCLUDES
   "${INCL_DIR}/SpidaCVX.h"
   "${INCL_DIR}/SpidaRVX.h"
   "${INCL_DIR}/SpidaCHEBX.h"
   "${INCL_DIR}/ChebInterpX.h"
   "${INCL_DIR}/SpidaCVT.h"
   "${INCL_DIR}/SpidaRVT.h"
   "${INCL_DIR}/SpidaR.h"
   "${INCL_DIR}/SpidaRCVT.h"
   "${INCL_DIR}/SpidaRRVT.h"
   "${INCL_DIR}/helper/constants.h"
   "${INCL_DIR}/helper/interp.h"
   "${INCL_DIR}/helper/funcs.hpp"
   "${INCL_DIR}/grid/grid.h"
   "${INCL_DIR}/grid/gridX.h"
   "${INCL_DIR}/grid/gridT.h"
   "${INCL_DIR}/grid/gridR.h"
   "${INCL_DIR}/grid/uniformX.h"
   "${INCL_DIR}/grid/uniformRVX.h"
   "${INCL_DIR}/grid/uniformCVX.h"
   "${INCL_DIR}/grid/uniformT.h"
   "${INCL_DIR}/grid/uniformRVT.h"
   "${INCL_DIR}/grid/uniformCVT.h"
   "${INCL_DIR}/grid/chebX.h"
   "${INCL_DIR}/grid/besselR.h"
   "${INCL_DIR}/propagator/propagator.h"
   "${INCL_DIR}/propagator/reporthandler.h"
   "${INCL_DIR}/shape/shape.h"
   "${INCL_DIR}/shape/shapeT.h"
   "${INCL_DIR}/shape/shapeR.h"
   "${INCL_DIR}/rkstiff/solver.h"
   "${INCL_DIR}/rkstiff/ETDAS.h"
   "${INCL_DIR}/rkstiff/ETDCS.h"
   "${INCL_DIR}/rkstiff/IFAS.h"
   "${INCL_DIR}/transform/chebX.h"
   "${INCL_DIR}/transform/fftCVT.h"
   "${INCL_DIR}/transform/fftRVT.h"
   "${INCL_DIR}/transform/fftCVX.h"
   "${INCL_DIR}/transform/fftRVX.h"
   "${INCL_DIR}/transform/hankelR.h"
   "${INCL_DIR}/transform/hankelfftRCVT.h"
   "${INCL_DIR}/transform/hankelfftRRVT.h"
)

set(BOOST_DIR "${SPIDA_SOURCE_DIR}/include/spida/boost_1_76_0/boost")
set(BOOST_INCLUDES 
	"${BOOST_DIR}/math/special_functions/airy.hpp"
	"${BOOST_DIR}/math/special_functions/bessel.hpp"
)

set(INCLUDE_FILES
	${BOOST_INCLUDES}
	${SPIDA_INCLUDES}
)

add_library(spida ${SRC_FILES} ${INCLUDE_FILES})
add_library(SPIDA::spida ALIAS spida)

target_include_directories(spida PUBLIC 
    $<BUILD_INTERFACE:${SPIDA_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
	$<BUILD_INTERFACE:${SPIDA_SOURCE_DIR}/include/spida/boost_1_76_0>
    # Need to find boost on local installs
	$<INSTALL_INTERFACE:include/spida/boost_1_76_0>
)

target_link_libraries(spida PUBLIC
    pw::pwutils
    kissfft
    nayukidct
)

set_target_properties(spida PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION 0
)

install(DIRECTORY ${INCL_DIR}
    DESTINATION include
)

#install(TARGETS spida pwutils nayukidct kissfft
install(TARGETS spida
	EXPORT SpidaTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
)

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

configure_file("${SPIDA_SOURCE_DIR}/cmake/SPIDAconfig.cmake.in"
    SpidaConfig.cmake @ONLY
)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfigVersion.cmake
    DESTINATION lib/cmake/spida
)







