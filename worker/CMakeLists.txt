# spida-worker: links SPIDA::spida directly, in-tree — no find_package/Conan
# sub-project round-trip. The worker was originally authored and verified as
# a separate Conan consumer of an installed spida package (see
# docs/adr/0002-worker-model-coverage-and-config-registry.md); moving it
# in-tree drops that whole layer, matching how test/ and demos/ already
# consume the library.
add_executable(spida-worker src/main.cpp)
target_link_libraries(spida-worker PRIVATE SPIDA::spida)
target_compile_features(spida-worker PRIVATE cxx_std_23)
