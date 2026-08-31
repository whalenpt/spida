
set(DEMO_TARGETS
    interp
    grid
    hankel
    transformT
    transformR
    transformRT
    transform_p_RT
    kdvCV
    kdvRV
    kdvPropagator
    ks
    ks_if34
    ks_if45dp
    burgers
    nlsR
    nlsRT
)

add_executable(interp interp_demo.cpp)
add_executable(grid grid_demo.cpp)
add_executable(hankel hankel.cpp)
add_executable(transformT transformT.cpp)
add_executable(transformR transformR.cpp)
add_executable(transformRT transformRT.cpp)
add_executable(transform_p_RT transform_p_RT.cpp)
add_executable(kdvCV kdvCV.cpp)
add_executable(kdvRV kdvRV.cpp)
add_executable(kdvPropagator kdv_propagator.cpp)
add_executable(ks ks.cpp)
add_executable(ks_if34 ks_if34.cpp)
add_executable(ks_if45dp ks_if45dp.cpp)
add_executable(burgers burgers.cpp)
add_executable(nlsR NLSR.cpp)
add_executable(nlsRT NLSRT.cpp)

foreach(target IN LISTS DEMO_TARGETS)
    target_link_libraries(${target} PUBLIC SPIDA::spida)
endforeach()



