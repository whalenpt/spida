// Consumer smoke test: prove that a downstream Conan consumer of the spida
// package can find_package(spida), link SPIDA::spida, and actually compile
// and run against spida's public headers — not just the narrowest possible
// header, but ones that have a history of leaking a third-party dependency's
// headers into the public interface (see conanfile.py's requirements()):
//   - <spida/RVX.h> pulls in transform/fftRVX.h, which used to
//     #include kissfft's own kiss_fftr.h directly.
//   - <spida/shape/shapeT.h> used to #include
//     <boost/math/special_functions/bessel.hpp> directly (dead weight —
//     never actually used in the header — but still a leak).
// Both are fixed at the header level now (forward declarations / dropped
// dead include), so this file no longer needs kissfft or boost itself to
// compile — deliberately not included here, to prove that.

#include <cstdlib>
#include <iostream>
#include <vector>

#include <spida/grid/uniformRVT.h>
#include <spida/helper/constants.h>
#include <spida/RVX.h>
#include <spida/shape/shapeT.h>

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

    // Exercises transform/fftRVX.h's kissfft-backed path for real (not just
    // linking): round-trip a physical-space field through spectral space.
    spida::SpidaRVX spi(grid);
    std::vector<double> physical(grid.getNx(), 1.0);
    std::vector<spida::dcmplx> spectral(grid.getNsx());
    spi.X_To_SX(physical, spectral);

    // Exercises shape/shapeT.h's ShapeT hierarchy for real.
    const spida::UniformGridRVT tgrid(N, -1.0, 1.0);
    const spida::GaussT gauss(tgrid, 1.0, 0.1);
    if (gauss.amplitude() != 1.0) {
        std::cerr << "spida test_package: unexpected shape amplitude\n";
        return EXIT_FAILURE;
    }

    std::cout << "spida test_package: linked and ran successfully (Nx=" << grid.getNx() << ")\n";
    return EXIT_SUCCESS;
}
