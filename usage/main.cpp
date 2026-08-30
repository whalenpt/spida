/*------------------------------------------------------------------------------
 *
 *    Description: Minimal example of consuming an installed SPIDA package —
 *    build a grid, run a transform round-trip, print the result. See
 *    README.md for the two-line CMake incantation this depends on
 *    (find_package(spida) + target_link_libraries(... SPIDA::spida)).
 *
------------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <vector>

#include <spida/grid/uniformRVX.h>
#include <spida/helper/constants.h>
#include <spida/RVX.h>

int main()
{
    constexpr unsigned N = 64;
    constexpr double a = 0.0;
    constexpr double b = 2.0 * spida::PI;

    const spida::UniformGridRVX grid(N, a, b);
    spida::SpidaRVX spi(grid);

    // A single cosine mode in physical space, round-tripped through spectral
    // space and back — cheap way to prove the library actually works, not
    // just that it links.
    const auto& x = grid.getX();
    std::vector<double> physical(grid.getNx());
    for (size_t i = 0; i < x.size(); i++)
        physical[i] = std::cos(x[i]);

    std::vector<spida::dcmplx> spectral(grid.getNsx());
    spi.X_To_SX(physical, spectral);

    std::vector<double> roundtrip(grid.getNx());
    spi.SX_To_X(spectral, roundtrip);

    std::cout << "SPIDA is installed and working (Nx=" << grid.getNx()
              << ", roundtrip[0]=" << roundtrip[0] << ")\n";

    return 0;
}
