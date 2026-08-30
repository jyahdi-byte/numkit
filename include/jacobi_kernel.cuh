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

__global__ void jacobi_kernel(double* d_old, double* d_new, CellType* d_types, double* d_faces_k, double* d_total_ks, unsigned char* d_active, int rows, int cols, int depth) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2) * (depth - 2)){return;}

    int plane = (rows - 2) * (cols - 2);
    int dep = index / plane;
    int rem = index % plane;
    int row = rem / (cols - 2);
    int col = rem % (cols - 2);

    int i = row + 1;
    int j = col + 1;
    int k = dep + 1;

    if (d_types[(k * rows * cols) + i * cols + j] == INTERIOR){
        d_new[(k * rows * cols) + i * cols + j] = update_cell(d_old, d_faces_k, d_total_ks, d_active, i, j, k, rows, cols, depth);
    }
}

#endif