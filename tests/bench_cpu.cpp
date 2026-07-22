#include <chrono>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "stats.hpp" 

void jacobi_sweep_n(Grid& g, int n){
    Grid newg = g;
    for (int k = 0; k < n; k++){
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                newg.at(i,j) = (g.at(i+1,j) + g.at(i-1,j) + g.at(i,j+1) + g.at(i,j-1))/4;
            }
        }
        g = newg;
    }
}

int main(){
    Grid g(100,100);
    for (int j = 0; j < g.getCols(); j++){
        g.at(0,j) = 100;
    }

    int rows = g.getRows();
    int cols = g.getCols();
    int sweeps = 10000;

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
}