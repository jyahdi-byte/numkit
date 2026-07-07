#include <iostream>
#include "grid.hpp"
#include "jacobi.hpp"
#include "gauss_seidel.hpp"

int main() {
    Grid g1(10, 10);
    for (int j = 0; j < g1.getCols(); j++) {
        g1.at(0, j) = 100.0;
    }
    Grid g2 = g1;   // identical copy, hot edge included

    int jac = jacobi_solve(g1, 1e-6, 10000);
    int gs  = gauss_seidel_solve(g2, 1e-6, 10000);

    std::cout << "Jacobi sweeps:       " << jac << "\n";
    std::cout << "Gauss-Seidel sweeps: " << gs  << "\n\n";

    std::cout << "            Jacobi        Gauss-Seidel\n";
    std::cout << "(1,5):  " << g1.at(1,5) << "      " << g2.at(1,5) << "\n";
    std::cout << "(5,5):  " << g1.at(5,5) << "      " << g2.at(5,5) << "\n";
    std::cout << "(8,5):  " << g1.at(8,5) << "      " << g2.at(8,5) << "\n";
    std::cout << "(5,2):  " << g1.at(5,2) << "      " << g2.at(5,2) << "\n";
    std::cout << "(5,7):  " << g1.at(5,7) << "      " << g2.at(5,7) << "\n";
}