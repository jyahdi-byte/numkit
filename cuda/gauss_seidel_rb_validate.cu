#include <iostream>
#include <cassert>
#include <cmath>

#include "grid.hpp"
#include "gauss_seidel_rb.hpp"
#include "gauss_seidel_rb_kernel.cuh"


int main(){
    Grid g(100,100);
    for (int j = 0; j < g.getCols(); j++){g.at(0,j) = 100;}
    g.maskRect(40, 40, 60, 60, HOLE);
    Grid g0 = g;;
 

    int rows = g.getRows();
    int cols = g.getCols();
    int sweeps = gauss_seidel_rb_solve(g, 1e-8, 50000);
    double* d;
    CellType* d_types;

    cudaMalloc((void**)&d, rows * cols * sizeof(double));
    cudaMalloc((void**)&d_types, rows * cols * sizeof(CellType));
    cudaMemcpy(d, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_types, g0.getTypesPtr(), rows * cols * sizeof(CellType), cudaMemcpyHostToDevice);


    int interior = (rows - 2) * (cols - 2); 
    int threadsPerBlock = 256;
    int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
    for (int i = 0; i < sweeps; i++){
        gs_rb_kernel_r<<<numBlocks, threadsPerBlock>>>(d, d_types, rows, cols);
        cudaDeviceSynchronize();
        gs_rb_kernel_b<<<numBlocks, threadsPerBlock>>>(d, d_types, rows, cols);
        cudaDeviceSynchronize();
    }
    double* g1 = (double*)malloc(rows * cols * sizeof(double));
    cudaMemcpy(g1, d, rows * cols * sizeof(double), cudaMemcpyDeviceToHost);

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            assert(std::abs(g.at(i,j) - g1[i * cols + j]) < 1e-9);
        }
    }

    return 0;
}