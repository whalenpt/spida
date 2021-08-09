cmake_minimum_required(VERSION 3.19)

set(THIS pwutils_test)
set(Sources
    mytest.cpp
)
add_executable(${THIS} ${Sources})
target_link_libraries(${THIS} PUBLIC
    gtest_main
    pwutils
)

add_test(
    NAME ${THIS}
    COMMAND ${THIS}
)
