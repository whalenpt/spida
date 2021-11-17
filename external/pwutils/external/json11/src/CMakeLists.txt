
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

install(DIRECTORY ${INCL_DIR}
	DESTINATION include
)

install(TARGETS json11
	EXPORT Json11Targets
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	RUNTIME DESTINATION bin
)

install(EXPORT Json11Targets
	FILE Json11Targets.cmake
	DESTINATION lib/cmake/json11
)

include(CMakePackageConfigHelpers)

configure_file(${json11_SOURCE_DIR}/cmake/Json11Config.cmake.in
    Json11Config.cmake @ONLY
)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/Json11Config.cmake
    DESTINATION lib/cmake/json11
)



