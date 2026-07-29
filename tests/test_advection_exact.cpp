#include <iostream>
#include <cassert>
#include <cmath>
#include "grid1d.hpp"
#include "advection.hpp"

int main(){

    Grid1D u(50, 0.02);
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}
    int iterations = 500;
    double cn = 1;

    SpaceTimeLog Log(iterations + 1, u.getPoints());
    advection_solve(u, iterations, 1, cn, Log);
    for (int i = 0; i < u.getPoints(); i++){
        int index = (int)(i - iterations * cn) % u.getPoints();
        if (index < 0){index = u.getPoints() + index;};
        if (index >= 5 && index <= 8) {assert(std::abs(u.at(i) - 1) < 1e-9);}
        else {assert(std::abs(u.at(i)) < 1e-9);}
    }

    std::cout << "PASS\n";

    return 0;
}