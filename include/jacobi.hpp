#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <cmath>
#include <vector>
#include <utility>

#include "grid.hpp"
#include "update_cell.hpp"

int jacobi_solve(Grid& g, double tol, int max_iter){
    std::vector<double> cur(g.getRows() * g.getCols()); std::vector<double> nxt(g.getRows() * g.getCols());
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            cur[i * g.getCols() + j] = g.at(i,j);
            nxt[i * g.getCols() + j] = g.at(i,j);
        }
    }

    for (int k = 0; k < max_iter; k++){
        double maxChange = 0;
        double change = 0;
        for (int i = 1; i < g.getRows() - 1; i++){ 
            for (int j = 1; j < g.getCols() - 1; j++){ 
                if (g.getType(i,j) == INTERIOR){
                    nxt[i * g.getCols() + j] = update_cell(cur.data(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                    change = std::abs(nxt[i * g.getCols() + j] - cur[i * g.getCols() + j]);
                    if (change >= maxChange){maxChange = change;}
                }
            }
        }
        std::swap(cur, nxt);
        if (maxChange <= tol){
            for (int i = 0; i < g.getRows(); i++){
                for (int j = 0; j < g.getCols(); j++){g.at(i,j) = cur[i * g.getCols() + j];}}
            return k+1;
        }
    }

     for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){g.at(i,j) = cur[i * g.getCols() + j];}}

    return max_iter;
}

#endif