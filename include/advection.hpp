#ifndef ADVECTION_HPP
#define ADVECTION_HPP

#include <iostream>
#include <cassert>
#include "grid1d.hpp"
#include "space_time_log.hpp"

double compute_dt(double cn, double dx, double c){
    assert(cn <= 1);
    return cn * dx / c;
}
 
double compute_cn(double dt, double dx, double c){ 
    return c * dt / dx;
}

int advection_solve(Grid1D& u, int num_steps, double c, double cn, SpaceTimeLog& log){
    if (cn > 1){
        std::cout << "cn = " << cn << ", exceeds limit of 1" << std::endl;
        assert(false);
    }
    Grid1D new_u(u.getPoints(), u.getDx());
    for (int k = 0; k < num_steps; k++){
        log.setRow(k, u);
        for (int i = 0; i < u.getPoints(); i++){
            new_u.at(i) = (1 - cn) * u.at(i) + cn * u.at(i-1);
        }
        u = new_u;
    }
    log.setRow(num_steps, u);
    return num_steps;
}

#endif