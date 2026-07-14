#include <iostream>
#include <cassert>
#include <cmath>
#include "grid.hpp"
#include "jacobi.hpp"

__global__ void jacobi_kernel(double* d_old, double* d_new, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2);
    int i = row + 1;
    int j = col + 1;

    d_new[i * cols + j] = (d_old[i * cols + j - 1] + d_old[i * cols + j + 1] + d_old[(i - 1) * cols + j] + d_old[(i + 1) * cols + j]) / 4;
}

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


    int interior = (rows - 2) * (cols - 2);
    int threadsPerBlock = 256;
    int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
    for (int i = 0; i < sweeps; i++){
        jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, rows, cols);
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