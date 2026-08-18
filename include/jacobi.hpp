#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <cmath>
#include <vector>

#include "grid.hpp"

int jacobi_solve(Grid& g, double tol, int max_iter){
    Grid newg = g;
    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                if (g.getType(i,j) == INTERIOR){
                    std::vector<double> jac_components = {g.at(i+1,j), g.at(i-1,j), g.at(i,j+1), g.at(i,j-1)}; 
                    if (g.getType(i+1,j) == HOLE){jac_components[0] = g.at(i,j);}
                    if (g.getType(i-1,j) == HOLE){jac_components[1] = g.at(i,j);}
                    if (g.getType(i,j+1) == HOLE){jac_components[2] = g.at(i,j);}
                    if (g.getType(i,j-1) == HOLE){jac_components[3] = g.at(i,j);}
                    newg.at(i,j) = (jac_components[0] + jac_components[1] + jac_components[2] + jac_components[3])/4;
                    change = std::abs(newg.at(i,j) - g.at(i,j));
                    if (change >= maxChange){
                        maxChange = change; 
                    }
                }
            }
        }
        g = newg;
        if (maxChange <= tol){
            return k+1;
        }
    }
    return max_iter;
}

#endif