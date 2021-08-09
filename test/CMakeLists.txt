cmake_minimum_required(VERSION 3.8)

set(THIS mytests)
set(SRC_FILES
	kiss_tests.cpp
	transform_tests.cpp
	dct_tests.cpp
	shape_tests.cpp
	nayukidct_tests.cpp
    fftw_tests.cpp
	fdct_tests.cpp
)

message("Run tests on library spida")

add_executable(${THIS} ${SRC_FILES})

message(STATUS "Link ${THIS} to pwutils and spida")
target_link_libraries(${THIS} PUBLIC
	gtest_main
    pw::pwutils
    SPIDA::spida
)

add_test(
	NAME ${THIS}
	COMMAND ${THIS}
)

