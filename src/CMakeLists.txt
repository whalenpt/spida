# -----------------------------
# Source files
# -----------------------------
set(SRC_FILES
    SpidaCVT.cpp
    SpidaRVT.cpp
    SpidaRCVT.cpp
    SpidaRRVT.cpp
    SpidaCVX.cpp
    SpidaRVX.cpp
    SpidaChebX.cpp
    SpidaR.cpp
    ChebInterpX.cpp
    helper/interp.cpp
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
	shape/shapeX.cpp
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

# -----------------------------
# Library target
# -----------------------------
add_library(spida ${SRC_FILES})
add_library(SPIDA::spida ALIAS spida)

# -----------------------------
# Includes
# -----------------------------
target_include_directories(spida PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# -----------------------------
# C++ standard
# -----------------------------
target_compile_features(spida PUBLIC cxx_std_17)

# -----------------------------
# Core dependencies
# -----------------------------
target_link_libraries(spida PUBLIC
    pw::pwutils
    kissfft::kissfft
    nayukidct
)

if(TARGET spdlog::spdlog)
    target_link_libraries(spida PUBLIC spdlog::spdlog)
endif()

# -----------------------------
# OpenBLAS optional
# -----------------------------
if(TARGET OpenBLAS::interface)
    target_link_libraries(spida PUBLIC OpenBLAS::interface)
endif()

# -----------------------------
# Boost handling
# -----------------------------
if(TARGET Boost::boost)
    target_link_libraries(spida PUBLIC Boost::boost)
else()
    message(STATUS "Using vendored Boost")

    set(LOCAL_BOOST_DIR "${PROJECT_SOURCE_DIR}/external/spidaboost_1_90_0")

    target_include_directories(spida PUBLIC
        $<BUILD_INTERFACE:${LOCAL_BOOST_DIR}>
    )

    install(DIRECTORY "${LOCAL_BOOST_DIR}/boost"
            DESTINATION include)
endif()

# -----------------------------
# Properties
# -----------------------------
set_target_properties(spida PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# -----------------------------
# Install
# -----------------------------
install(TARGETS spida
    EXPORT SpidaTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION include
)

install(EXPORT SpidaTargets
    FILE SpidaTargets.cmake
    NAMESPACE SPIDA::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spida
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/SPIDAconfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfig.cmake
    @ONLY
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SpidaConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spida
)
