#ifndef BLACK_SCHOLES_HPP
#define BLACK_SCHOLES_HPP

#include <cmath>
#include <cassert>
#include <iostream>

#include "grid1d.hpp"
#include "space_time_log.hpp"
#include "diffusion.hpp"

int black_scholes_solve(Grid1D& g, int num_steps, double x_max, double rate, double dtau, double sigma, double K, double alpha, double dx){
    double r = compute_r(std::pow(sigma, 2) / 2, dtau, dx);
    if (r < 0 || r > 0.5){
        if (r < 0){std::cout << "r = "<< r << ", below minimum of 0" << std::endl;}
        else{std::cout << "r = " << r << ", exceeds maximum of 0.5" << std::endl;}
        assert(false);
    }
    Grid1D newg = g;
    for (int k = 0; k < num_steps; k++){
        newg.at(0) = 0;
        newg.at(g.getPoints() - 1) = (std::exp(x_max + (dtau * k) * rate) - K);
        for (int i = 1; i < g.getPoints() - 1; i++){
            newg.at(i) = g.at(i) + r * (g.at(i+1) - 2 * g.at(i) + g.at(i-1)) + (rate - alpha) * dtau * (g.at(i) - g.at(i-1)) / dx;
        }
        g = newg; 
    }
    return num_steps;
}

#endif
