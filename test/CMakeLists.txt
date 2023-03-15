cmake_minimum_required(VERSION 3.15)


set(GRID_FILES grid_tests.cpp)
set(SHAPE_FILES shape_tests.cpp)
set(TRANSFORM_FILES dct_tests.cpp 
	            fftT_tests.cpp
                fftX_tests.cpp
                hankel_tests.cpp
                hankelfftRRVT_tests.cpp
                hankelfftRCVT_tests.cpp)

set(SPIDA_FILES ${GRID_FILES} ${SHAPE_FILES} ${TRANSFORM_FILES})

set(KISS_FILES kiss_tests.cpp)
set(NAYUKI_FILES nayuki_tests.cpp)
	            

message("Run tests on grids")
add_executable(gridtest ${GRID_FILES})
message("Run tests on shapes")
add_executable(shapetest ${SHAPE_FILES})
message("Run tests on transforms")
add_executable(transformtest ${TRANSFORM_FILES})

message("Run tests on library spida")
add_executable(spidatest ${SPIDA_FILES})

message("Run tests on kissfft")
add_executable(kisstest ${KISS_FILES})
message("Run tests on Nayuki FDCT")
add_executable(nayukitest ${NAYUKI_FILES})

message(STATUS "Link spidatest to pwutils and spida")
target_link_libraries(shapetest PUBLIC
    gtest_main
    SPIDA::spida
)

target_link_libraries(gridtest PUBLIC
    gtest_main
    SPIDA::spida
)

target_link_libraries(transformtest PUBLIC
    gtest_main
    SPIDA::spida
    pw::pwutils
)

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
	NAME gridtest
	COMMAND gridtest
)
add_test(
	NAME shapetest
	COMMAND shapetest
)
add_test(
	NAME transformtest
	COMMAND transformtest
)
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
