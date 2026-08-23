#include <iostream>

#include "grid.hpp"
#include "conjugate_gradient.hpp"

int main(){
    std::cout << "=== Single hole ===\n";
    for (int N = 10; N <= 80; N *= 2){
        if (N < 20) continue;
        Grid g1(N, N);
        for (int j = 0; j < g1.getCols(); j++) {g1.at(0, j) = 10.0;}
        g1.maskRect(N/3, N/3, N/3 + 3, N/3 + 3, HOLE);
        Grid g2 = g1;

        int cg = cg_solve(g1, 1e-10, 10000);
        int pcg = pcg_solve(g2, 1e-10, 10000);
        std::cout << "N = " << N << ": CG = " << cg << ", PCG = " << pcg;
        std::cout << " (" << (double)cg / pcg << "x)\n";
    }

    std::cout << "\n=== Checkerboard of holes ===\n";
    for (int N = 10; N <= 80; N *= 2){
        Grid g1(N, N);
        for (int j = 0; j < g1.getCols(); j++) {g1.at(0, j) = 10.0;}
        for (int i = 1; i < N - 1; i += 4){
            for (int j = 1; j < N - 1; j += 4){
                g1.maskRect(i, j, i + 3, j + 3, HOLE);
            }
        }
        Grid g2 = g1;

        int cg = cg_solve(g1, 1e-10, 10000);
        int pcg = pcg_solve(g2, 1e-10, 10000);
        std::cout << "N = " << N << ": CG = " << cg << ", PCG = " << pcg;
        std::cout << " (" << (double)cg / pcg << "x)\n";
    }

    return 0;
}