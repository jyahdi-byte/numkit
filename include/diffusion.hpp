#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP

#include <cassert>
#include <iostream>
#include <cmath>

#include "grid1d.hpp"
#include "space_time_log.hpp"

double compute_r(double alpha, double dt, double dx){
    return (alpha * dt / std::pow(dx, 2));
}

int diffusion_solve(Grid1D& u, int num_steps, double r, SpaceTimeLog& log){
    if (r < 0 || r > 0.5){
        if (r < 0){std::cout << "r = "<< r << ", below minimum of 0" << std::endl;}
        else{std::cout << "r = " << r << ", exceeds maximum of 0.5" << std::endl;}
        assert(false);
    }
    Grid1D new_u(u.getPoints(), u.getDx());
    for (int k = 0; k < num_steps; k++){
        log.setRow(k, u);
        for (int i = 0; i < u.getPoints(); i++){
            new_u.at(i) = u.at(i) + r * (u.at(i+1) - 2 * u.at(i) + u.at(i-1));
        }
        u = new_u;
    }
    log.setRow(num_steps, u);
    return num_steps;
}


#endif