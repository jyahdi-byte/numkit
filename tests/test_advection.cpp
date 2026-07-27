#include <iostream>
#include <cassert>
#include <cmath>
#include "grid1d.hpp"
#include "advection.hpp"

int main(){
    Grid1D u(20, 0.02);
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}
    double sum_before = u.sum();
    u.print(); std::cout << "Sum: " << sum_before << "\n\n";

    advection_solve(u, 30, 1, 0.01);
    double sum_after = u.sum();
    u.print(); std::cout << "Sum: " << sum_after << "\n\n";

    assert(std::abs(sum_after - sum_before) < 1e-9);
    std::cout << "PASS\n";

    return 0;
}