// Minimal consumer smoke test: prove that a downstream Conan consumer of the
// spida package can find_package(spida), link SPIDA::spida, and actually run
// something that touches spida code (not just link it).

#include <cstdlib>
#include <iostream>

#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>

int main()
{
    constexpr unsigned N = 64;
    constexpr double a = 0.0;
    constexpr double b = 2.0 * spida::PI;

    const spida::UniformGridRVX grid(N, a, b);

    if (grid.getNx() != N) {
        std::cerr << "spida test_package: unexpected grid size\n";
        return EXIT_FAILURE;
    }

    std::cout << "spida test_package: linked and ran successfully (Nx=" << grid.getNx() << ")\n";
    return EXIT_SUCCESS;
}
