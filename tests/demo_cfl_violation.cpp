#include <cassert>
#include "grid1d.hpp"
#include "advection.hpp"

int main(){
    //This is a Deliberate CFL Violation

    Grid1D u(50, 0.02);
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}

    int iterations = 500;
    SpaceTimeLog Log(iterations + 1, u.getPoints());
    advection_solve(u, iterations, 1, 100, Log);

    return 0;
} 