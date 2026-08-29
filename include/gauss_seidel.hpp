#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP

#include <cmath>
#include <vector>

#include "grid.hpp"
#include "update_cell.hpp"

int gauss_seidel_solve(Grid& g, double tol, int max_iter){
    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
               if (g.getType(i,j) == INTERIOR){
                    double oldPoint = g.at(i,j); 
                    g.at(i,j) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                    change = std::abs(g.at(i,j) - oldPoint);
                    if (change >= maxChange){
                        maxChange = change; 
                    }
                }
            }
        }
        if (maxChange <= tol){
            return k+1;
        }
    }
    return max_iter;
}

int gauss_seidel_solve(Grid3D& g, double tol, int max_iter){
    for (int n = 0; n < max_iter; n++){
        double maxChange = 0;
        double change = 0;
        for (int k = 1; k < g.getDepth() - 1; k++){
            for (int i = 1; i < g.getRows() - 1; i++){
                for (int j = 1; j < g.getCols() - 1; j++){
                if (g.getType(i,j,k) == INTERIOR){
                        double oldPoint = g.at(i,j,k); 
                        g.at(i,j,k) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, k, g.getRows(), g.getCols(), g.getDepth());
                        change = std::abs(g.at(i,j,k) - oldPoint);
                        if (change >= maxChange){
                            maxChange = change; 
                        }
                    }
                }
            }
        }
        if (maxChange <= tol){
            return n+1;
        }
    }
    return max_iter;
}

#endif