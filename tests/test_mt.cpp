#include <iostream>
#include <chrono>
#include "grid.hpp"
#include "jacobi_mt.hpp"
#include "jacobi.hpp"

int main() {
    Grid g1(200, 200);
    Grid g2(200, 200);

    for (int j = 0; j < g1.getCols(); j++) {
        g1.at(0, j) = 100.0;
    }
    for (int j = 0; j < g2.getCols(); j++) {
        g2.at(0, j) = 100.0;
    }

    // Timing scaffold (chrono) AI-assisted; solver code and test design are mine.
    auto t0 = std::chrono::steady_clock::now();
    int sweeps1 = jacobi_mt_solve(g1, 1e-6, 100000, 4);
    auto t1 = std::chrono::steady_clock::now();
    int sweeps2 = jacobi_solve(g2, 1e-6, 100000);
    auto t2 = std::chrono::steady_clock::now();

    double mt_time = std::chrono::duration<double>(t1 - t0).count();
    double st_time = std::chrono::duration<double>(t2 - t1).count();

    std::cout << "With Threads (4): \n";
    std::cout << "Converged in " << sweeps1 << " sweeps, " << mt_time << " s\n";
    std::cout << "Near hot edge  (20,100): " << g1.at(20, 100) << "\n";
    std::cout << "Center         (100,100): " << g1.at(100, 100) << "\n";
    std::cout << "Far side       (160,100): " << g1.at(160, 100) << "\n";
    std::cout << "Left of center (100,80): " << g1.at(100, 80) << "\n";
    std::cout << "Right of center(100,140): " << g1.at(100, 140) << "\n";

    std::cout << "\nWithout Threads: \n";
    std::cout << "Converged in " << sweeps2 << " sweeps, " << st_time << " s\n";
    std::cout << "Near hot edge  (20,100): " << g2.at(20, 100) << "\n";
    std::cout << "Center         (100,100): " << g2.at(100, 100) << "\n";
    std::cout << "Far side       (160,100): " << g2.at(160, 100) << "\n";
    std::cout << "Left of center (100,80): " << g2.at(100, 80) << "\n";
    std::cout << "Right of center(100,140): " << g2.at(100, 140) << "\n";

    std::cout << "\nSpeedup: " << st_time / mt_time << "x\n";
    return 0;
}