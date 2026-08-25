#ifndef SOR_RB_HPP
#define SOR_RB_HPP


#include <cmath>
#include <numbers>
#include <vector>

#include "grid.hpp"
#include "update_cell.hpp"


int sor_rb_solve(Grid& g, double tol, int max_iter, double omega){
    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 0){
                    double oldPoint = g.at(i,j);
                    double avg = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                    double delta = omega * (avg - oldPoint);
                    g.at(i,j) = oldPoint + delta;
                    change = std::abs(delta);
                    if (change >= maxChange){
                        maxChange = change;
                    }
                }
            }
        }
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 1){
                    double oldPoint = g.at(i,j);
                    double avg = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                    double delta = omega * (avg - oldPoint);
                    g.at(i,j) = oldPoint + delta;
                    change = std::abs(delta);
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

int sor_rb_solve(Grid& g, double tol, int max_iter) {
    double h = 1.0 / (g.getRows() - 1);
    double omega = 2.0 / (1.0 + std::sin(std::numbers::pi * h));
    return sor_rb_solve(g, tol, max_iter, omega);
}

#endif