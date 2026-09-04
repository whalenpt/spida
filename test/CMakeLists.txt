# ------------------------------------------------------------------------------
# Test source groups
# ------------------------------------------------------------------------------

set(GRID_FILES
    grid_tests.cpp
    grid_extended_tests.cpp
)

set(SHAPE_FILES
    shape_tests.cpp
)

set(TRANSFORM_FILES
    dct_tests.cpp
    fftT_tests.cpp
    fftX_tests.cpp
    hankel_tests.cpp
    hankelfftRRVT_tests.cpp
    hankelfftRCVT_tests.cpp
    spidaR_tests.cpp
    chebInterp_tests.cpp
    spida_wrappers_tests.cpp
)

set(SPIDA_FILES
    ${GRID_FILES}
    ${SHAPE_FILES}
    ${TRANSFORM_FILES}
)

set(KISS_FILES kiss_tests.cpp)
set(NAYUKI_FILES nayuki_tests.cpp)
set(INTERP_FILES interp_tests.cpp)
set(SOLVER_FILES solver_tests.cpp)
set(PROPAGATOR_FILES propagator_tests.cpp)
set(CONFIG_FILES config_tests.cpp)

# ------------------------------------------------------------------------------
# Helper function to reduce repetition
# ------------------------------------------------------------------------------

function(add_spida_test target)
    add_executable(${target} ${ARGN})

    target_link_libraries(${target} PUBLIC
        GTest::gtest_main
        SPIDA::spida
    )

    add_test(NAME ${target} COMMAND ${target})
endfunction()

# ------------------------------------------------------------------------------
# Core SPIDA tests
# ------------------------------------------------------------------------------

add_spida_test(gridtest ${GRID_FILES})
add_spida_test(shapetest ${SHAPE_FILES})

add_spida_test(transformtest ${TRANSFORM_FILES})
target_link_libraries(transformtest PRIVATE nayukidct)

add_spida_test(spidatest ${SPIDA_FILES})
target_link_libraries(spidatest PRIVATE nayukidct)

# ------------------------------------------------------------------------------
# External libraries
# ------------------------------------------------------------------------------

add_executable(kisstest ${KISS_FILES})
target_link_libraries(kisstest PUBLIC
    GTest::gtest_main
    kissfft::kissfft
)
target_include_directories(kisstest PRIVATE ${PROJECT_SOURCE_DIR}/src)
add_test(NAME kisstest COMMAND kisstest)

add_executable(nayukitest ${NAYUKI_FILES})
target_link_libraries(nayukitest PUBLIC
    GTest::gtest_main
    nayukidct
)
target_include_directories(nayukitest PRIVATE ${PROJECT_SOURCE_DIR}/src)
add_test(NAME nayukitest COMMAND nayukitest)

# ------------------------------------------------------------------------------
# Additional SPIDA components
# ------------------------------------------------------------------------------

add_spida_test(interptest ${INTERP_FILES})

add_spida_test(solvertest ${SOLVER_FILES})
add_spida_test(propagatortest ${PROPAGATOR_FILES})
add_spida_test(utilstest utils_tests.cpp)
add_spida_test(configtest ${CONFIG_FILES})
