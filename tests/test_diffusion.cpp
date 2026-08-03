#include <iostream>
#include <cassert>
#include <cmath>

#include "grid1d.hpp"
#include "space_time_log.hpp"
#include "diffusion.hpp"
#include "ppm.hpp"
 
int main(){  
    Grid1D u(50, 0.02);
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}
    for (int i = 27; i < 42; i++){u.at(i) = 0.8;}
    double sum_before = u.sum();
    u.print(); std::cout << "Sum: " << sum_before << "\n\n";

    int iterations = 500;
    SpaceTimeLog Log(iterations + 1, u.getPoints());
    diffusion_solve(u, iterations, 0.33, Log);
    double sum_after = u.sum();
    u.print(); std::cout << "Sum: " << sum_after << "\n\n";

    assert(std::abs(sum_after - sum_before) < 1e-9);
    std::cout << "PASS\n";

    write_ppm(Log, "DiffusionMap.ppm"); 
 
    return 0;
} 