#include <iostream>
#include <cassert>
#include <cmath>

#include "grid.hpp"
#include "jacobi.hpp"
#include "gauss_seidel.hpp"
#include "sor.hpp"
#include "jacobi_mt.hpp"
#include "gauss_seidel_rb.hpp"
#include "sor_rb.hpp"
#include "conjugate_gradient.hpp"

int main() {
    Grid3D g1(10, 10, 10);
    for (int i = 0; i < g1.getRows(); i++){
        for (int j = 0; j < g1.getCols(); j++){
            for (int k = 0; k < g1.getDepth(); k++){
                if (i == 0 || j == 0 || k == 0 || i == g1.getRows() - 1 || j == g1.getCols() - 1 || k == g1.getDepth() - 1){
                    g1.at(i,j,k) = i;
                }
            }
        }
    }
    g1.maskRect(4,4,4,6,6,6,HOLE);
    Grid3D g2 = g1;
    Grid3D g3 = g1;
    Grid3D g4 = g1;
    Grid3D g5 = g1;
    Grid3D g6 = g1;
    Grid3D g7 = g1;
    Grid3D g8 = g1;
    Grid3D g0 = g1;

    int jac = jacobi_solve(g1, 1e-10, 10000);
    int gs  = gauss_seidel_solve(g2, 1e-10, 10000);
    int sor = sor_solve(g3, 1e-10, 10000, 1.5);
    int sor1 = sor_solve(g4, 1e-10, 10000);
    int gs_rb = gauss_seidel_rb_solve(g5, 1e-10, 10000);
    int sor_rb = sor_rb_solve(g6, 1e-10, 10000);
    int cg = cg_solve(g7, 1e-10, 10000);
    int pcg = pcg_solve(g8, 1e-10, 10000);

    std::cout << "Jacobi sweeps:                   " << jac << "\n";
    std::cout << "Gauss-Seidel sweeps:             " << gs  << "\n";
    std::cout << "SOR (w=1.5) sweeps:              " << sor << "\n";
    std::cout << "Sor (Theoretical Best w) sweeps: " << sor1 << "\n";
    std::cout << "Gauss-Seidel-RB sweeps:                " << gs_rb << "\n";
    std::cout << "Sor-RB (Theoretical Best w) sweeps: " << sor_rb << "\n";
    std::cout << "Conjugate Gradient sweeps:       " << cg << "\n";
    std::cout << "Preconditioned Conjugate Gradient sweeps: " << pcg << "\n\n";

    std::cout << "(3,6,6)  Jacobi: " << g1.at(3,6,6) << "\n"
              << "   GS: "        << g2.at(3,6,6) << "\n"
              << "   SOR(w=1.5): "       << g3.at(3,6,6) << "\n"
              << "   SOR(Theoretical Best W): " << g4.at(3,6,6) << "\n"
              << "   GS-RB: " << g5.at(3,6,6) << "\n"
              << "   SOR-RB: " << g6.at(3,6,6) << "\n"
              << "   CG: " << g7.at(3,6,6) << "\n"
              << "   PCG: " << g8.at(3,6,6) << "\n\n";

    assert(std::abs(g1.at(3,6,6) - g2.at(3,6,6)) < 1e-6);
    assert(std::abs(g2.at(3,6,6) - g3.at(3,6,6)) < 1e-6);
    assert(std::abs(g3.at(3,6,6) - g4.at(3,6,6)) < 1e-6);
    assert(std::abs(g4.at(3,6,6) - g5.at(3,6,6)) < 1e-6);
    assert(std::abs(g5.at(3,6,6) - g6.at(3,6,6)) < 1e-6);
    assert(std::abs(g6.at(3,6,6) - g7.at(3,6,6)) < 1e-6);
    assert(std::abs(g7.at(3,6,6) - g8.at(3,6,6)) < 1e-6);
    assert(g1.at(5,5,5) == g2.at(5,5,5) && g2.at(5,5,5) == g3.at(5,5,5) && g3.at(5,5,5) == g4.at(5,5,5) && g4.at(5,5,5) == g5.at(5,5,5)
            && g5.at(5,5,5) == g6.at(5,5,5) && g6.at(5,5,5) == g7.at(5,5,5) && g7.at(5,5,5) == g8.at(5,5,5) && g8.at(5,5,5) == g0.at(5,5,5));

    std::cout << "PASS\n";
}