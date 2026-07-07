#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP
#include "grid.hpp"
#include <cmath>

int gauss_seidel_solve(Grid& g, double tol, int max_iter){
    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double oldPoint = g.at(i,j);
                g.at(i,j) = (g.at(i+1,j) + g.at(i-1,j) + g.at(i,j+1) + g.at(i,j-1))/4;
                change = std::abs(g.at(i,j) - oldPoint);
                if (change >= maxChange){
                    maxChange = change;
                }
            }
        }
        if (maxChange <= tol){
            return k+1;
        }
    }
    return max_iter;
}

#endif