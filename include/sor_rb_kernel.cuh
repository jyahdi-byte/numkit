#ifndef SOR_RB_KERNEL_CUH
#define SOR_RB_KERNEL_CUH

#include "grid.hpp"

__global__ void sor_rb_kernel_r(double* d, CellType* d_types, int rows, int cols, double omega) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1; 

    if (d_types[i * cols + j] == INTERIOR && (i + j) % 2 == 0){
        double oldPoint = d[i * cols + j];
        double components[4] = {d[i * cols + j - 1], d[i * cols + j + 1], d[(i - 1) * cols + j], d[(i + 1) * cols + j]};
        if (d_types[i * cols + j - 1] == HOLE){components[0] = d[i * cols + j];}
        if (d_types[i * cols + j + 1] == HOLE){components[1] = d[i * cols + j];}
        if (d_types[(i - 1) * cols + j] == HOLE){components[2] = d[i * cols + j];}
        if (d_types[(i + 1) * cols + j] == HOLE){components[3] = d[i * cols + j];}
        double avg = (components[0] + components[1] + components[2] + components[3])/4;
        double delta = omega * (avg - oldPoint);
        d[i * cols + j] = oldPoint + delta;
    }
}

__global__ void sor_rb_kernel_b(double* d, CellType* d_types, int rows, int cols, double omega) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1;

    if (d_types[i * cols + j] == INTERIOR && (i + j) % 2 == 1){
        double oldPoint = d[i * cols + j];
        double components[4] = {d[i * cols + j - 1], d[i * cols + j + 1], d[(i - 1) * cols + j], d[(i + 1) * cols + j]};
        if (d_types[i * cols + j - 1] == HOLE){components[0] = d[i * cols + j];}
        if (d_types[i * cols + j + 1] == HOLE){components[1] = d[i * cols + j];}
        if (d_types[(i - 1) * cols + j] == HOLE){components[2] = d[i * cols + j];}
        if (d_types[(i + 1) * cols + j] == HOLE){components[3] = d[i * cols + j];}
        double avg = (components[0] + components[1] + components[2] + components[3])/4;
        double delta = omega * (avg - oldPoint);
        d[i * cols + j] = oldPoint + delta;
    }
}
 
#endif