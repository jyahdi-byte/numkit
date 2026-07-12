#include <iostream>
#include "grid.hpp"
#include "sor.hpp"

int main(){
    Grid g1(100, 100);
    for (int j = 0; j < g1.getCols(); j++){
        g1.at(0, j) = 100.0;
    }
    Grid g2 = g1;

    int sweeps_auto = sor_solve(g1, 1e-6, 100000);
    int sweeps_manual = sor_solve(g2, 1e-6, 100000, 1.5);

    std::cout << "Auto omega sweeps:   " << sweeps_auto << "\n";
    std::cout << "Manual omega (1.5):  " << sweeps_manual << "\n";
    std::cout << "(1,50) auto:   " << g1.at(1, 50) << "\n";
    std::cout << "(1,50) manual: " << g2.at(1, 50) << "\n";

    return 0;
}