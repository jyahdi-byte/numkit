#ifndef JACOBI_KERNEL_CUH
#define JACOBI_KERNEL_CUH

#include "grid.hpp"
#include "update_cell.hpp"  

__global__ void jacobi_kernel(double* d_old, double* d_new, CellType* d_types, double* d_faces_k, double* d_total_ks, unsigned char* d_active, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1;

    if (d_types[i * cols + j] == INTERIOR){
        d_new[i * cols + j] = update_cell(d_old, d_faces_k, d_total_ks, d_active, i, j, rows, cols);
    }
}

#endif