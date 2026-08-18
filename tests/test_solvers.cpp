#include <iostream>
#include <cassert>
#include <cmath>

#include "grid.hpp"
#include "jacobi.hpp"
#include "gauss_seidel.hpp"
#include "sor.hpp"
#include "jacobi_mt.hpp"

int main() {
    Grid g1(10, 10);
    for (int j = 0; j < g1.getCols(); j++) {g1.at(0, j) = 100.0;}
    g1.maskRect(4,4,6,6,HOLE);
    Grid g2 = g1;
    Grid g3 = g1;
    Grid g4 = g1;
    Grid g5 = g1;
    Grid g0 = g1;

    int jac = jacobi_solve(g1, 1e-10, 10000);
    int gs  = gauss_seidel_solve(g2, 1e-10, 10000);
    int sor = sor_solve(g3, 1e-10, 10000, 1.5);
    int sor1 = sor_solve(g4, 1e-10, 10000);
    int jac_mt = jacobi_mt_solve(g5, 1e-10, 10000, 4);

    std::cout << "Jacobi sweeps:                   " << jac << "\n";
    std::cout << "Gauss-Seidel sweeps:             " << gs  << "\n";
    std::cout << "SOR (w=1.5) sweeps:              " << sor << "\n";
    std::cout << "Sor (Theoretical Best w) sweeps: " << sor1 << "\n";
    std::cout << "Jacobi_mt sweeps:                " << jac_mt << "\n";

    std::cout << "(3,6)  Jacobi: " << g1.at(3,6)
              << "   GS: "        << g2.at(3,6)
              << "   SOR(w=1.5): "       << g3.at(3,6) 
              << "   SOR(Theoretical Best W): " << g4.at(3,6) 
              << "   Jacobi_mt: " << g5.at(3,6) << "\n";

    assert(std::abs(g1.at(3,6) - g2.at(3,6)) < 1e-6);
    assert(std::abs(g2.at(3,6) - g3.at(3,6)) < 1e-6);
    assert(std::abs(g3.at(3,6) - g4.at(3,6)) < 1e-6);
    assert(std::abs(g4.at(3,6) - g5.at(3,6)) < 1e-6);
    assert(g1.at(5,5) == g2.at(5,5) && g2.at(5,5) == g3.at(5,5) && g3.at(5,5) == g4.at(5,5) && g4.at(5,5) == g5.at(5,5)
            && g5.at(5,5) == g0.at(5,5));

    std::cout << "PASS\n";
}