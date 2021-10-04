
cmake_minimum_required(VERSION 3.19)

add_executable(interp_ex interp_example.cpp)
add_executable(grid_ex grid_example.cpp)
add_executable(hankel_ex hankel_example.cpp)
add_executable(transformT_ex transformT_example.cpp)
add_executable(transformR_ex transformR_example.cpp)
add_executable(transformRT_ex transformRT_example.cpp)
add_executable(kdvCV_ex kdv_example.cpp)
add_executable(kdvRV_ex kdv_exampleRV.cpp)
add_executable(kdvPropagator_ex kdv_example_propagator.cpp)

target_link_libraries(interp_ex PUBLIC
    SPIDA::spida
    pw::pwutils
    ParamBin::parambin
)

target_link_libraries(grid_ex PUBLIC SPIDA::spida)
target_link_libraries(hankel_ex PUBLIC SPIDA::spida)

target_link_libraries(transformRT_ex PUBLIC 
    SPIDA::spida
    pw::pwutils
)
target_link_libraries(transformT_ex PUBLIC 
    SPIDA::spida
    pw::pwutils
)
target_link_libraries(transformR_ex PUBLIC 
    SPIDA::spida
    pw::pwutils
)

target_link_libraries(kdvPropagator_ex PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(kdvCV_ex PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(kdvRV_ex PUBLIC SPIDA::spida pw::pwutils)



