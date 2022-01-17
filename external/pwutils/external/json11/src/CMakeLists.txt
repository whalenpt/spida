

cmake_minimum_required(VERSION 3.19)

set(INCL_DIR "${json11_SOURCE_DIR}/include/json11")
add_library(json11 json11.cpp ${INCL_DIR}/json11.hpp)
set_target_properties(json11 PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(NOT WIN32)
    target_compile_options(json11
      PRIVATE -fno-rtti -fno-exceptions -Wall)
endif()

target_include_directories(json11 PRIVATE
	$<BUILD_INTERFACE:${json11_SOURCE_DIR}/include/json11>
	$<INSTALL_INTERFACE:include/json11>
)

target_include_directories(json11 PUBLIC 
	$<BUILD_INTERFACE:${json11_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

# Set warning flags, which may vary per platform

install(TARGETS json11 
    EXPORT JSON11Targets
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	RUNTIME DESTINATION bin
)
install(DIRECTORY ${INCL_DIR}
	DESTINATION include
)
install(EXPORT JSON11Targets
    FILE JSON11Targets.cmake
    DESTINATION lib/cmake/json11
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${json11_SOURCE_DIR}/cmake/JSON11Config.cmake.in
    ${json11_BINARY_DIR}/cmake/JSON11Config.cmake
    INSTALL_DESTINATION lib/cmake/json11
    PATH_VARS CMAKE_INSTALL_LIBDIR
)

install(
    FILES "${json11_BINARY_DIR}/cmake/JSON11Config.cmake"
    DESTINATION lib/cmake/json11
)


