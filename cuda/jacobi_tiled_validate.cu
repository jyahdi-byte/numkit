#include <cassert>
#include <cmath>
#include "grid.hpp"
#include "jacobi.hpp"
#include "jacobi_tiled_kernel.cuh"

int main(){
    Grid g(100,100);
    Grid g0(100,100);
    for (int j = 0; j < g.getCols(); j++){
        g.at(0,j) = 100;
        g0.at(0,j) = 100;
    }

    int rows = g.getRows();
    int cols = g.getCols();
    int sweeps = jacobi_solve(g, 1e-8, 50000);
    double* d_old;
    double* d_new;
    cudaMalloc((void**)&d_old, rows * cols * sizeof(double));
    cudaMalloc((void**)&d_new, rows * cols * sizeof(double));
    cudaMemcpy(d_old, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_new, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16,16);
    dim3 numBlocks((rows + 16 - 1) / 16, (cols + 16 - 1) / 16);

    for (int i = 0; i < sweeps; i++){
        jacobi_tiled_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, rows, cols);
        cudaDeviceSynchronize();

        double* temp = d_old;
        d_old = d_new;
        d_new = temp;
    }
    double* g1 = (double*)malloc(rows * cols * sizeof(double));
    cudaMemcpy(g1, d_old, rows * cols * sizeof(double), cudaMemcpyDeviceToHost);

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            assert(std::abs(g.at(i,j) - g1[i * cols + j]) < 1e-9);
        }
    }

    return 0;
}