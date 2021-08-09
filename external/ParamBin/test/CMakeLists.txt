cmake_minimum_required(VERSION 3.19)

set(Sources
    parambin_test.cpp
)

add_executable(parambin_test ${Sources})
target_link_libraries(parambin_test PUBLIC
    gtest_main
    parambin
)

add_test(
    NAME parambin_test
    COMMAND parambin_test
)
