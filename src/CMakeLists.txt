# -----------------------------
# Source files
# -----------------------------
set(SRC_FILES
    utils/math.cpp
    utils/stats.cpp
    utils/strings.cpp
    utils/threads.cpp
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
target_include_directories(spida
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

# -----------------------------
# C++ standard
# -----------------------------
target_compile_features(spida PUBLIC cxx_std_23)

# -----------------------------
# Core dependencies
# -----------------------------
target_link_libraries(spida PUBLIC
    kissfft::kissfft
    Threads::Threads
)

target_link_libraries(spida PRIVATE
    nayukidct
)

if(TARGET spdlog::spdlog)
    target_link_libraries(spida PRIVATE spdlog::spdlog)
endif()

# -----------------------------
# nlohmann_json
# -----------------------------
find_package(nlohmann_json CONFIG QUIET)
if(nlohmann_json_FOUND)
    target_link_libraries(spida PUBLIC nlohmann_json::nlohmann_json)
else()
    message(STATUS "nlohmann_json not found via Conan; fetching single header")
    set(_json_header "${CMAKE_CURRENT_BINARY_DIR}/nlohmann/json.hpp")
    if(NOT EXISTS "${_json_header}")
        file(DOWNLOAD
            "https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp"
            "${_json_header}"
            SHOW_PROGRESS
        )
    endif()
    target_include_directories(spida PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
    )
    install(FILES "${_json_header}" DESTINATION include/nlohmann)
endif()

# -----------------------------
# OpenBLAS optional
# -----------------------------
if(TARGET OpenBLAS::interface)
    target_link_libraries(spida PRIVATE OpenBLAS::interface)
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

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/utils
    DESTINATION include
)

install(EXPORT SpidaTargets
    FILE SpidaTargets.cmake
    NAMESPACE SPIDA::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spida
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/spidaConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/SPIDAconfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/spidaConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spida
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/spidaConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/spidaConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spida
)
