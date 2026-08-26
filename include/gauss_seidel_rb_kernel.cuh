#ifndef GAUSS_SEIDEL_RB_KERNEL_CUH
#define GAUSS_SEIDEL_RB_KERNEL_CUH

#include "grid.hpp"
#include "update_cell.hpp"

__global__ void gs_rb_kernel_r(double* d, CellType* d_types, double* d_faces_k, double* d_total_ks, unsigned char* d_active, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1; 

    if (d_types[i * cols + j] == INTERIOR && (i + j) % 2 == 0){
        d[i * cols + j] = update_cell(d, d_faces_k, d_total_ks, d_active, i, j, rows, cols);
    }
}

__global__ void gs_rb_kernel_b(double* d, CellType* d_types, double* d_faces_k, double* d_total_ks, unsigned char* d_active, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);   
    int col = index % (cols - 2); 
    int i = row + 1;
    int j = col + 1;

    if (d_types[i * cols + j] == INTERIOR && (i + j) % 2 == 1){
        d[i * cols + j] = update_cell(d, d_faces_k, d_total_ks, d_active, i, j, rows, cols);
    }
}
 
#endif