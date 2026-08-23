#ifndef GAUSS_SEIDEL_RB_HPP
#define GAUSS_SEIDEL_RB_HPP
#include "grid.hpp"
#include <cmath>
#include <vector>
 
int gauss_seidel_rb_solve(Grid& g, double tol, int max_iter){
    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
               if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 0){
                    double oldPoint = g.at(i,j);
                    std::vector<double> components = {g.at(i+1,j), g.at(i-1,j), g.at(i,j+1), g.at(i,j-1)}; 
                    if (g.getType(i+1,j) == HOLE){components[0] = g.at(i,j);}
                    if (g.getType(i-1,j) == HOLE){components[1] = g.at(i,j);}
                    if (g.getType(i,j+1) == HOLE){components[2] = g.at(i,j);}
                    if (g.getType(i,j-1) == HOLE){components[3] = g.at(i,j);}
                    g.at(i,j) = (components[0] + components[1] + components[2] + components[3])/4;
                    change = std::abs(g.at(i,j) - oldPoint);
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
                    std::vector<double> components = {g.at(i+1,j), g.at(i-1,j), g.at(i,j+1), g.at(i,j-1)}; 
                    if (g.getType(i+1,j) == HOLE){components[0] = g.at(i,j);}
                    if (g.getType(i-1,j) == HOLE){components[1] = g.at(i,j);}
                    if (g.getType(i,j+1) == HOLE){components[2] = g.at(i,j);}
                    if (g.getType(i,j-1) == HOLE){components[3] = g.at(i,j);}
                    g.at(i,j) = (components[0] + components[1] + components[2] + components[3])/4;
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

#endif