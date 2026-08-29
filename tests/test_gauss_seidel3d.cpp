 #include <iostream>
#include <cassert>
#include <cmath>

#include "grid.hpp"
#include "gauss_seidel.hpp"

int main() {
    Grid3D g(10, 10, 10); 

    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            for (int k = 0; k < g.getDepth(); k++){
                if (i == 0 || i == g.getRows() - 1 || j == 0 || j == g.getCols() - 1 || k == 0 || k == g.getDepth() - 1){
                    g.at(i,j,k) = i;
                }
            }
        }
    }

    int sweeps = gauss_seidel_solve(g, 1e-8, 10000);
    std::cout << "Converged in " << sweeps << " sweeps\n";

    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            for (int k = 0; k < g.getDepth(); k++){
                assert(std::abs(i - g.at(i,j,k)) < 1e-6);
            }
        }
    }

    std::cout << "PASS\n";
    return 0;
}