#include <iostream>
#include <string>
#include <cmath>

#include "grid.hpp"
#include "sor.hpp"
#include "ppm.hpp"

int main(int argc, char* argv[]){
    Grid g(200, 300);
    for (int j = 0; j < g.getCols(); j++){
        g.at(0,j) = 100.0;
    }

    bool obstacle = false;
    for (int a = 1; a < argc; a++){
        if (std::string(argv[a]) == "--obstacle"){
            obstacle = true;
        }
    }

    if (obstacle){
        // circular insulated obstacle, centered in the grid
        g.maskCircle(g.getRows() / 2, g.getCols() / 2, 40, HOLE);
    }

    int sweeps = sor_solve(g, 1e-6, 50000, 1.5);
    std::cout << sweeps << "\n";
    write_ppm_composite(g, "heat.ppm", 0.5);

    return 0;
}