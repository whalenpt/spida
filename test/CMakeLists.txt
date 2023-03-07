cmake_minimum_required(VERSION 3.8)

set(SPIDA_FILES
	dct_tests.cpp
	fftX_tests.cpp
	fftT_tests.cpp
	hankel_tests.cpp
    hankelfftRRVT_tests.cpp
    hankelfftRCVT_tests.cpp
	shape_tests.cpp
	grid_tests.cpp
)

set(KISS_FILES kiss_tests.cpp)
set(NAYUKI_FILES nayuki_tests.cpp)

message("Run tests on library spida")
add_executable(spidatest ${SPIDA_FILES})
message("Run tests on kissfft")
add_executable(kisstest ${KISS_FILES})
message("Run tests on Nayuki FDCT")
add_executable(nayukitest ${NAYUKI_FILES})

message(STATUS "Link spida test to pwutils and spida")
target_link_libraries(spidatest PUBLIC
	gtest_main
    SPIDA::spida
    pw::pwutils
)

target_link_libraries(kisstest PUBLIC
	gtest_main
    kissfft
    pw::pwutils
)

target_link_libraries(nayukitest PUBLIC
	gtest_main
	nayukidct
    pw::pwutils
)

#set(BOOST_DIR "${Spida_SOURCE_DIR}/external/boost_1_76_0")
#target_include_directories(${THIS} PRIVATE ${BOOST_DIR})

add_test(
	NAME spidatest
	COMMAND spidatest
)
add_test(
	NAME kisstest
	COMMAND kisstest
)
add_test(
	NAME nayukitest
	COMMAND nayukitest
)
