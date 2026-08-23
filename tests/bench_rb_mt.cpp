#include <iostream>
#include <chrono>
#include <cmath>

#include "grid.hpp"
#include "gauss_seidel_rb.hpp"
#include "gauss_seidel_rb_mt.hpp"
#include "sor_rb.hpp"
#include "sor_rb_mt.hpp"

int main() {

    Grid warmup(200, 200);
        for (int j = 0; j < warmup.getCols(); j++) {
            warmup.at(0, j) = 100.0;
        }
        gauss_seidel_rb_solve(warmup, 1e-6, 100000);

    std::cout << "=== Gauss-Seidel Red-Black ===\n";

    Grid gs_single(200, 200);
    for (int j = 0; j < gs_single.getCols(); j++) {
        gs_single.at(0, j) = 100.0;
    }
    // Timing scaffold (chrono) AI-assisted; solver code and test design are mine.
    auto s0 = std::chrono::steady_clock::now();
    int single_sweeps = gauss_seidel_rb_solve(gs_single, 1e-6, 100000);
    auto s1 = std::chrono::steady_clock::now();
    double single_time = std::chrono::duration<double>(s1 - s0).count();
    std::cout << "Single-threaded: converged in " << single_sweeps << " sweeps, " << single_time << " s\n";

    for (int i = 1; i < 9; i *= 2){
        Grid g1(200, 200);
        for (int j = 0; j < g1.getCols(); j++) {
            g1.at(0, j) = 100.0;
        }

        auto t0 = std::chrono::steady_clock::now();
        int mt_sweeps = gauss_seidel_rb_mt_solve(g1, 1e-6, 100000, i);
        auto t1 = std::chrono::steady_clock::now();
        double mt_time = std::chrono::duration<double>(t1 - t0).count();

        std::cout << "With Threads (" << i << ")  : \n";
        std::cout << "Converged in " << mt_sweeps << " sweeps, " << mt_time << " s";
        std::cout << ", speedup vs single-threaded: " << (single_time / mt_time) << "x\n";
    }

    std::cout << "\n=== SOR Red-Black ===\n";

    Grid sor_single(200, 200);
    for (int j = 0; j < sor_single.getCols(); j++) {
        sor_single.at(0, j) = 100.0;
    }
    auto r0 = std::chrono::steady_clock::now();
    int sor_single_sweeps = sor_rb_solve(sor_single, 1e-6, 100000);
    auto r1 = std::chrono::steady_clock::now();
    double sor_single_time = std::chrono::duration<double>(r1 - r0).count();
    std::cout << "Single-threaded: converged in " << sor_single_sweeps << " sweeps, " << sor_single_time << " s\n";

    for (int i = 1; i < 9; i *= 2){
        Grid g1(200, 200);
        for (int j = 0; j < g1.getCols(); j++) {
            g1.at(0, j) = 100.0;
        }

        auto t0 = std::chrono::steady_clock::now();
        int mt_sweeps = sor_rb_mt_solve(g1, 1e-6, 100000, i);
        auto t1 = std::chrono::steady_clock::now();
        double mt_time = std::chrono::duration<double>(t1 - t0).count();

        std::cout << "With Threads (" << i << ")  : \n";
        std::cout << "Converged in " << mt_sweeps << " sweeps, " << mt_time << " s";
        std::cout << ", speedup vs single-threaded: " << (sor_single_time / mt_time) << "x\n";
    }

    return 0;
}