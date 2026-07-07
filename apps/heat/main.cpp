#include <iostream>
#include "grid.hpp"
#include "sor.hpp"
#include "ppm.hpp"

int main(){
    Grid g(200, 300);
    for (int j = 0; j < g.getCols(); j++){
        g.at(0,j) = 100.0;
    }

    int sweeps = sor_solve(g, 1e-6, 50000, 1.5);
    std::cout << sweeps << "\n";
    write_ppm(g, "heat.ppm");

    return 0;
}