#ifndef SOR_HPP
#define SOR_HPP
#include "grid.hpp"
#include <cmath>

double PI = 3.14159265358979;

int sor_solve(Grid& g, double tol, int max_iter, double omega){
    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double oldPoint = g.at(i,j);
                double avg = (g.at(i+1,j) + g.at(i-1,j) + g.at(i,j+1) + g.at(i,j-1))/4;
                double delta = omega * (avg - oldPoint);
                g.at(i,j) = oldPoint + delta;
                change = std::abs(delta);
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

int sor_solve(Grid& g, double tol, int max_iter) {
    double h = 1.0 / (g.getRows() - 1);
    double omega = 2.0 / (1.0 + std::sin(PI * h));
    return sor_solve(g, tol, max_iter, omega);
}

#endif