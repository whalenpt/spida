
cmake_minimum_required(VERSION 3.19)

add_executable(interp interp_demo.cpp)
add_executable(grid grid_demo.cpp)
add_executable(hankel hankel.cpp)
add_executable(transformT transformT.cpp)
add_executable(transformR transformR.cpp)
add_executable(transformRT transformRT.cpp)
add_executable(kdvCV kdvCV.cpp)
add_executable(kdvRV kdvRV.cpp)
add_executable(kdvPropagator kdv_propagator.cpp)
add_executable(ks ks.cpp)
add_executable(burgers burgers.cpp)
add_executable(nlsR NLSR.cpp)

target_link_libraries(interp PUBLIC
    SPIDA::spida
    pw::pwutils
    ParamBin::parambin
)

target_link_libraries(grid PUBLIC SPIDA::spida)
target_link_libraries(hankel PUBLIC SPIDA::spida)
target_link_libraries(transformRT PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(transformT PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(transformR PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(kdvPropagator PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(kdvCV PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(kdvRV PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(ks PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(burgers PUBLIC SPIDA::spida pw::pwutils)
target_link_libraries(nlsR PUBLIC SPIDA::spida pw::pwutils)



