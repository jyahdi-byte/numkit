#include <cassert>
#include <cmath>
#include "grid.hpp"
#include "jacobi.hpp"
#include "jacobi_tiled_kernel.cuh" 
#include "cuda_check.cuh"

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
    double* d_faces_k; 
    double* d_total_ks;
    unsigned char* d_active;
    CUDA_CHECK(cudaMalloc((void**)&d_old, rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_new, rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 4 * rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_total_ks, rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_active, rows * cols * sizeof(unsigned char)));
    CUDA_CHECK(cudaMemcpy(d_old, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_new, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_faces_k, g0.getFacesKPtr(), 4 * rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_total_ks, g0.getTotalKPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_active, g0.getActivePtr(), rows * cols * sizeof(unsigned char), cudaMemcpyHostToDevice));

    dim3 threadsPerBlock(16,16);
    dim3 numBlocks((rows + threadsPerBlock.x - 1) / threadsPerBlock.x, (cols + threadsPerBlock.x - 1) / threadsPerBlock.x);

    for (int i = 0; i < sweeps; i++){
        jacobi_tiled_kernel<<<numBlocks, threadsPerBlock, (threadsPerBlock.x + 2) * (threadsPerBlock.x + 2) * sizeof(double)>>>(d_old, d_new, d_faces_k, d_total_ks, d_active, rows, cols);
        CUDA_CHECK(cudaDeviceSynchronize());

        double* temp = d_old;
        d_old = d_new;
        d_new = temp;
    }
    double* g1 = (double*)malloc(rows * cols * sizeof(double));
    CUDA_CHECK(cudaMemcpy(g1, d_old, rows * cols * sizeof(double), cudaMemcpyDeviceToHost));

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            assert(std::abs(g.at(i,j) - g1[i * cols + j]) < 1e-9);
        }
    }

    free(g1);
    CUDA_CHECK(cudaFree(d_old));
    CUDA_CHECK(cudaFree(d_new));
    CUDA_CHECK(cudaFree(d_faces_k));
    CUDA_CHECK(cudaFree(d_total_ks));
    CUDA_CHECK(cudaFree(d_active));

    return 0;
}