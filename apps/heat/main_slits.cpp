#include <iostream>
#include <cmath>

#include "grid.hpp"
#include "sor.hpp"
#include "ppm.hpp"

int main(){
    Grid g(80, 300);
    for (int j = 0; j < g.getCols(); j++){
        g.at(0,j) = 100.0;
    }

    int wallRowStart = 30;
    int wallRowEnd   = 34;   // wall is 4 rows thick

    int slit1Start = 90,  slit1End = 100;
    int slit2Start = 200, slit2End = 210;

    g.maskRect(wallRowStart, 0,            wallRowEnd, slit1Start, HOLE);
    g.maskRect(wallRowStart, slit1End,     wallRowEnd, slit2Start, HOLE);
    g.maskRect(wallRowStart, slit2End,     wallRowEnd, g.getCols(), HOLE);

    int sweeps = sor_solve(g, 1e-6, 50000, 1.5);
    std::cout << sweeps << "\n";
    write_ppm_composite(g, "heat_slits.ppm", 0.5);

    return 0;
}