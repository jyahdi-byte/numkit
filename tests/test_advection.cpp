#include <iostream>
#include <cassert>
#include <cmath>
#include "grid1d.hpp"
#include "advection.hpp"
#include "ppm.hpp"

int main(){ 
    Grid1D u(50, 0.02);
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}
    double sum_before = u.sum();
    u.print(); std::cout << "Sum: " << sum_before << "\n\n";

    int iterations = 500;
    SpaceTimeLog Log(iterations + 1, u.getPoints());
    advection_solve(u, iterations, 1, 0.5, Log);
    double sum_after = u.sum();
    u.print(); std::cout << "Sum: " << sum_after << "\n\n";

    assert(std::abs(sum_after - sum_before) < 1e-9);
    std::cout << "PASS\n";

    write_ppm(Log, "AdvectionMap.ppm"); 
 
    return 0;
} 