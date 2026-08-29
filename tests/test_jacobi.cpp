#include <iostream>
#include <cassert>
#include <cmath>

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

    assert(g.at(1,5) >= g.at(5,5));
    assert(g.at(5,5) >= g.at(8,5));
    assert(std::abs(g.at(5,2) - g.at(5,7)) < 1e-9);
    assert(std::abs(g.at(5,4) - g.at(5,5)) < 1e-9); 

    Grid g1(10,10);
    for (int j = 0; j < g1.getCols(); j++){
        g1.at(0,j) = 5.0;
    }

    g1.maskRect(4,4,6,6, HOLE);
    Grid oldg1 = g1;
    jacobi_solve(g1, 1e-6, 1000);
    std::cout << "diff near hole: " << std::abs(g1.at(3,4) - g1.at(2,4)) << "\n";
    for (int j = 0; j < 5; j++){assert(std::abs(g1.at(1,j) - g1.at(1,9-j)) < 1e-9);}
    for (int i = 4; i < 6; i++){for (int j = 4; j < 6; j++){assert(g1.at(i,j) == oldg1.at(i,j));}}
    std::cout << "diff near hole: " << std::abs(g1.at(3,4) - g1.at(2,4)) << "\n";
    assert(std::abs(g1.at(3,4) - g1.at(2,4)) < 1.0);

 

    std::cout << "PASS\n";
    return 0;
}