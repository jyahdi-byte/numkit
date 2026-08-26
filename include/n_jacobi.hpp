#ifndef N_JACOBI_HPP
#define N_JACOBI_HPP

#include <vector>
#include <utility>

#include "grid.hpp"
#include "update_cell.hpp"

void jacobi_sweep_n(Grid& g, int n){
    int rows = g.getRows();
    int cols = g.getCols();
    std::vector<double> cur(rows * cols);
    std::vector<double> nxt(rows * cols);
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            cur[i * cols + j] = g.at(i, j);
            nxt[i * cols + j] = g.at(i, j);
        }
    }

    for (int k = 0; k < n; k++){
        for (int i = 1; i < rows - 1; i++){
            for (int j = 1; j < cols - 1; j++){
                if (g.getType(i, j) == INTERIOR){
                    nxt[i * cols + j] = update_cell(cur.data(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, rows, cols);
                }
            }
        }
        std::swap(cur, nxt);
    }

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            g.at(i, j) = cur[i * cols + j];
        }
    }
}

#endif