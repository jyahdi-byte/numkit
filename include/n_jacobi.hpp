#ifndef N_JACOBI_HPP
#define N_JACOBI_HPP

#include "grid.hpp"

void jacobi_sweep_n(Grid& g, int n){
    Grid newg = g;
    for (int k = 0; k < n; k++){
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                newg.at(i,j) = (g.at(i+1,j) + g.at(i-1,j) + g.at(i,j+1) + g.at(i,j-1))/4;
            }
        }
        g = newg;
    }
}

#endif