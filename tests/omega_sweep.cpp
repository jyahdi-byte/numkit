#include <iostream>
#include "grid.hpp"
#include "sor.hpp"

int main(){
    double omega = 0.95;
    for (int i = 0; i < 19; i++){
        Grid g1(10, 10);
        for (int j = 0; j < g1.getCols(); j++) {
            g1.at(0, j) = 100.0;
        }
        omega += 0.05;
        int sweeps = sor_solve(g1, 1e-6, 10000, omega); 
        std::cout << "omega = " << omega << "  sweeps = " << sweeps << "\n";
    }
}

