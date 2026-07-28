#ifndef ADVECTION_HPP
#define ADVECTION_HPP

#include <cassert>
#include "grid1d.hpp"
#include "space_time_log.hpp"

int advection_solve(Grid1D& u, int num_steps, double c, double dt, SpaceTimeLog& log){
    double cn = c * dt / u.getDx();
    assert(cn <= 1);
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