#include <iostream>

#include "grid1d.hpp"
#include "space_time_log.hpp"
#include "diffusion.hpp"
 
int main(){  
    //This is a deliberate r violation

    Grid1D u(50, 0.02);
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}
    for (int i = 27; i < 42; i++){u.at(i) = 0.8;}

    int iterations = 500;
    SpaceTimeLog Log(iterations + 1, u.getPoints());
    diffusion_solve(u, iterations, 5, Log);
 
    return 0;
} 