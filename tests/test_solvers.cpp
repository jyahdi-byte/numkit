#include <iostream>
#include "grid.hpp"
#include "jacobi.hpp"
#include "gauss_seidel.hpp"
#include "sor.hpp"

int main() {
    Grid g1(10, 10);
    for (int j = 0; j < g1.getCols(); j++) {
        g1.at(0, j) = 100.0;
    }
    Grid g2 = g1;
    Grid g3 = g1;

    int jac = jacobi_solve(g1, 1e-6, 10000);
    int gs  = gauss_seidel_solve(g2, 1e-6, 10000);
    int sor = sor_solve(g3, 1e-6, 10000, 1.5);

    std::cout << "Jacobi sweeps:       " << jac << "\n";
    std::cout << "Gauss-Seidel sweeps: " << gs  << "\n";
    std::cout << "SOR (w=1.5) sweeps:  " << sor << "\n\n";

    std::cout << "(5,5)  Jacobi: " << g1.at(5,5)
              << "   GS: "        << g2.at(5,5)
              << "   SOR: "       << g3.at(5,5) << "\n";
}