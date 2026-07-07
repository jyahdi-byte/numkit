#include "jacobi.hpp"

#include <iostream>
#include "grid.hpp"
#include "jacobi.hpp"

int main() {
    Grid g(10, 10);

    for (int j = 0; j < g.getCols(); j++) {
        g.at(0, j) = 100.0;   // top edge held hot
    }

    int sweeps = jacobi_solve(g, 1e-6, 10000);

    std::cout << "Converged in " << sweeps << " sweeps\n";
    std::cout << "Near hot edge  (1,5): " << g.at(1, 5) << "\n";
    std::cout << "Center         (5,5): " << g.at(5, 5) << "\n";
    std::cout << "Far side       (8,5): " << g.at(8, 5) << "\n";
    std::cout << "Left of center (5,2): " << g.at(5, 2) << "\n";
    std::cout << "Right of center(5,7): " << g.at(5, 7) << "\n";
    return 0;
}