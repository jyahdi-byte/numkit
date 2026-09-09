#include <chrono>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "stats.hpp" 
#include "n_jacobi.hpp"

int main(){
    Grid3D g(40,40,40);
    for (int j = 0; j < g.getCols(); j++){
        for (int k = 0; k < g.getDepth(); k++){
            g.at(0,j,k) = 100;
        }
    }

    int rows = g.getRows();
    int cols = g.getCols();
    int depth = g.getDepth();
    int sweeps = 5000;

    jacobi_sweep_n(g, 10);
    std::vector<double> times;

    for (int k = 0; k < 10; k++){

        // AI added: timer around the real timed loop
        auto t1 = std::chrono::high_resolution_clock::now();

            jacobi_sweep_n(g, sweeps);

        // AI added: stop timer, convert to ms
        auto t2 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        std::cout << "Elapsed: " << ms << " ms for " << sweeps << " sweeps\n";
        times.push_back(ms);
    }
    std::cout << "\n";
    std::cout << "Mean Elapsed Time: " << mean(times) << "\n";
    std::cout << "Elapsed Time Standard Deviation: " << dev(times) << "\n";
    
    return 0;
}