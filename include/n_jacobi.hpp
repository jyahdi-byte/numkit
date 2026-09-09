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

void jacobi_sweep_n(Grid3D& g, int n){
    int rows = g.getRows();
    int cols = g.getCols();
    int depth = g.getDepth();
    std::vector<double> cur(rows * cols * depth);
    std::vector<double> nxt(rows * cols * depth);
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            for (int k = 0; k < depth; k++){
                cur[i * cols * depth + j * depth + k] = g.at(i, j, k);
                nxt[i * cols * depth + j * depth + k] = g.at(i, j, k);
            }
        }
    }

    for (int l = 0; l < n; l++){
        for (int i = 1; i < rows - 1; i++){
            for (int j = 1; j < cols - 1; j++){
                for (int k = 1; k < depth - 1; k++){
                    if (g.getType(i, j, k) == INTERIOR){
                        nxt[i * cols * depth + j * depth + k] = update_cell(cur.data(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, k, rows, cols, depth);
                    }
                }
            }
        }
        std::swap(cur, nxt);
    }

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            for (int k = 0; k < depth; k++){
                g.at(i, j, k) = cur[i * cols * depth + j * depth + k];
            }
        }
    }
}

#endif