cmake_minimum_required(VERSION 3.19)

set(THIS mytests)
set(SRC_FILES
    test.cpp
    string_test.cpp
    report_test.cpp
    read_test.cpp
)

message("Run tests on library pwutils")

add_executable(${THIS} ${SRC_FILES})
target_link_libraries(${THIS} PUBLIC
    gtest_main
    pwutils
)

add_test(
    NAME ${THIS}
    COMMAND ${THIS}
)

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

