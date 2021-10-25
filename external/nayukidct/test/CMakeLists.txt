cmake_minimum_required(VERSION 3.19)
set(CMAKE_CXX_STANDARD 11)

add_executable(nayukitest FastDctTest.cpp)
target_link_libraries(nayukitest PUBLIC nayukidct)

add_executable(dcttest dcttest.cpp)
target_link_libraries(dcttest PUBLIC 
    nayukidct
    gtest_main
)

add_test(
	NAME dcttest
	COMMAND dcttest
)


