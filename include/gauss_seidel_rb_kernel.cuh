#ifndef GAUSS_SEIDEL_RB_KERNEL_CUH
#define GAUSS_SEIDEL_RB_KERNEL_CUH

#include "grid.hpp"

__global__ void gs_rb_kernel_r(double* d, CellType* d_types, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1;

    if (d_types[i * cols + j] == INTERIOR && (i + j) % 2 == 0){
        double components[4] = {d[i * cols + j - 1], d[i * cols + j + 1], d[(i - 1) * cols + j], d[(i + 1) * cols + j]};
        if (d_types[i * cols + j - 1] == HOLE){components[0] = d[i * cols + j];}
        if (d_types[i * cols + j + 1] == HOLE){components[1] = d[i * cols + j];}
        if (d_types[(i - 1) * cols + j] == HOLE){components[2] = d[i * cols + j];}
        if (d_types[(i + 1) * cols + j] == HOLE){components[3] = d[i * cols + j];}
        d[i * cols + j] = (components[0] + components[1] + components[2] + components[3]) / 4;
    }
}

__global__ void gs_rb_kernel_b(double* d, CellType* d_types, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1;

    if (d_types[i * cols + j] == INTERIOR && (i + j) % 2 == 1){
        double components[4] = {d[i * cols + j - 1], d[i * cols + j + 1], d[(i - 1) * cols + j], d[(i + 1) * cols + j]};
        if (d_types[i * cols + j - 1] == HOLE){components[0] = d[i * cols + j];}
        if (d_types[i * cols + j + 1] == HOLE){components[1] = d[i * cols + j];}
        if (d_types[(i - 1) * cols + j] == HOLE){components[2] = d[i * cols + j];}
        if (d_types[(i + 1) * cols + j] == HOLE){components[3] = d[i * cols + j];}
        d[i * cols + j] = (components[0] + components[1] + components[2] + components[3]) / 4;
    }
}
 
#endif