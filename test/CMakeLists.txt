cmake_minimum_required(VERSION 3.8)

set(THIS mytests)
set(SRC_FILES
	kiss_tests.cpp
	nayukidct_tests.cpp
	dct_tests.cpp
	fft_tests.cpp
	transform_tests.cpp
	shape_tests.cpp
	grid_tests.cpp
)

#set(BOOST_INCLUDES 
#	"${BOOST_DIR}/boost/math/special_functions/airy.hpp"
#	"${BOOST_DIR}/boost/math/special_functions/bessel.hpp"
#)

message("Run tests on library spida")

add_executable(${THIS} ${SRC_FILES})

message(STATUS "Link ${THIS} to pwutils and spida")
target_link_libraries(${THIS} PUBLIC
	gtest_main
    pw::pwutils
    SPIDA::spida
    fftw3
)

set(BOOST_DIR "${Spida_SOURCE_DIR}/external/boost_1_76_0")
target_include_directories(${THIS} PRIVATE ${BOOST_DIR})

add_test(
	NAME ${THIS}
	COMMAND ${THIS}
)

