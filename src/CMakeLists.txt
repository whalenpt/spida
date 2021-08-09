

cmake_minimum_required(VERSION 3.19.1)

set(SRC_FILES
    grid/chebX.cpp
    grid/uniformX.cpp
	grid/uniformT.cpp
	interp.cpp
	model/modelT.cpp
	propagator/propagatorRT_RC.cpp
	propagator/propagatorRTp_RC.cpp
	propagator/propagatorT.cpp
	report/reportcenter.cpp
	report/reporthandler.cpp
	shape/shapeT.cpp
	solver/solver.cpp
	solver/solverETDAS.cpp
	solver/solverETDCS.cpp
	spidaPeriodicT.cpp
    spidaChebX.cpp
    spidaFFTWX.cpp
    spidaInterpFFTWX.cpp
	transform/periodicT.cpp
	transform/chebFFTWX.cpp
	transform/chebX.cpp
    #	transform/fftwT.cpp
    #	transform/transformX.cpp
)

set(INCL_DIR "${Spida_SOURCE_DIR}/include/spida")
set(SPIDA_INCLUDES
   "${INCL_DIR}/constants.h"
   "${INCL_DIR}/grid/uniformT.h"
   "${INCL_DIR}/grid/gridX.h"
   "${INCL_DIR}/grid/chebX.h"
   "${INCL_DIR}/grid/uniformX.h"
   "${INCL_DIR}/interp.h"
   "${INCL_DIR}/model/model.h"
   "${INCL_DIR}/model/modelT.h"
   "${INCL_DIR}/propagator/propagator.h"
   "${INCL_DIR}/propagator/propagatorT.h"
   "${INCL_DIR}/report/reportcenter.h"
   "${INCL_DIR}/report/reporthandler.h"
   "${INCL_DIR}/shape/shape.hpp"
   "${INCL_DIR}/shape/shapeT.h"
   "${INCL_DIR}/solver/solver.h"
   "${INCL_DIR}/solver/solverETDAS.h"
   "${INCL_DIR}/solver/solverETDCS.h"
   "${INCL_DIR}/spidaT.h"
   "${INCL_DIR}/spidaX.h"
   "${INCL_DIR}/spidaChebX.h"
   "${INCL_DIR}/spidaFFTWX.h"
   "${INCL_DIR}/spidaInterpFFTWX.h"
   "${INCL_DIR}/transform/transformT.h"
   "${INCL_DIR}/transform/transformX.h"
   "${INCL_DIR}/transform/periodicT.h"
   "${INCL_DIR}/transform/chebX.h"
   "${INCL_DIR}/transform/chebFFTWX.h"

   #   "${INCL_DIR}/transform/fftwT.h"
   #"${INCL_DIR}/transform/transformX.h"
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
    PUBLIC fftw3
)

set_target_properties(spida PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION 0
)

#install(TARGETS nlohmann_json spida
install(TARGETS spida
	EXPORT SpidaTargets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
)

#install(FILES ${SPIDA_INCLUDES}
#	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/spida"
#)

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







