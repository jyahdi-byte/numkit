#include <iostream>
#include <chrono>
#include "grid.hpp"
#include "jacobi_mt.hpp"
#include "jacobi.hpp"

int main() {
    for (int i = 1; i < 9; i *= 2){

        Grid g1(200, 200);
        for (int j = 0; j < g1.getCols(); j++) {
            g1.at(0, j) = 100.0;
        }

        // Timing scaffold (chrono) AI-assisted; solver code and test design are mine.
        auto t0 = std::chrono::steady_clock::now();
        int sweeps1 = jacobi_mt_solve(g1, 1e-6, 100000, i);
        auto t1 = std::chrono::steady_clock::now();


        double mt_time = std::chrono::duration<double>(t1 - t0).count();

        std::cout << "With Threads (" << i << ")  : \n";
        std::cout << "Converged in " << sweeps1 << " sweeps, " << mt_time << " s\n";
        std::cout << "Near hot edge  (20,100): " << g1.at(20, 100) << "\n";
        std::cout << "Center         (100,100): " << g1.at(100, 100) << "\n";
        std::cout << "Far side       (160,100): " << g1.at(160, 100) << "\n";
        std::cout << "Left of center (100,80): " << g1.at(100, 80) << "\n";
        std::cout << "Right of center(100,140): " << g1.at(100, 140) << "\n";

    }

    Grid g1(200,200);
    Grid g2(200,200);
    auto c0 = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000; i++){
        g1 = g2;
    }
    auto c1 = std::chrono::steady_clock::now();
    double perCopy = std::chrono::duration<double>(c1 - c0).count() / 1000.0;
    std::cout << "\nPer copy: " << perCopy << " s\n";
    std::cout << "Copy total over 68427 sweeps: " << perCopy * 68427 << " s\n";

    return 0;
}