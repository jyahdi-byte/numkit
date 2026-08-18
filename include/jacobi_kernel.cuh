#ifndef JACOBI_KERNEL_CUH
#define JACOBI_KERNEL_CUH

#include "grid.hpp"

__global__ void jacobi_kernel(double* d_old, double* d_new, CellType* d_types, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1;

    if (d_types[i * cols + j] == INTERIOR){
        double components[4] = {d_old[i * cols + j - 1], d_old[i * cols + j + 1], d_old[(i - 1) * cols + j], d_old[(i + 1) * cols + j]};
        if (d_types[i * cols + j - 1] == HOLE){components[0] = d_old[i * cols + j];}
        if (d_types[i * cols + j + 1] == HOLE){components[1] = d_old[i * cols + j];}
        if (d_types[(i - 1) * cols + j] == HOLE){components[2] = d_old[i * cols + j];}
        if (d_types[(i + 1) * cols + j] == HOLE){components[3] = d_old[i * cols + j];}
        d_new[i * cols + j] = (components[0] + components[1] + components[2] + components[3]) / 4;
    }
}

#endif