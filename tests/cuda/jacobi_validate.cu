#include <iostream>
#include <cassert>
#include <cmath>
#include "grid.hpp"
#include "jacobi.hpp"
#include "jacobi_kernel.cuh"


int main(){
    Grid g(100,100);
    for (int j = 0; j < g.getCols(); j++){g.at(0,j) = 100;}
    g.maskRect(40, 40, 60, 60, HOLE);
    Grid g0 = g;;


    int rows = g.getRows();
    int cols = g.getCols();
    int sweeps = jacobi_solve(g, 1e-8, 50000);
    double* d_old;
    double* d_new;
    CellType* d_types;
    double* d_faces_k; 
    double* d_total_ks; 
    unsigned char* d_active;
    cudaMalloc((void**)&d_old, rows * cols * sizeof(double));
    cudaMalloc((void**)&d_new, rows * cols * sizeof(double));
    cudaMalloc((void**)&d_types, rows * cols * sizeof(CellType));
    cudaMalloc((void**)&d_faces_k, 4 * rows * cols * sizeof(double));
    cudaMalloc((void**)&d_total_ks, rows * cols * sizeof(double));
    cudaMalloc((void**)&d_active, rows * cols * sizeof(unsigned char));
    cudaMemcpy(d_old, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_new, g0.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_types, g0.getTypesPtr(), rows * cols * sizeof(CellType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_faces_k, g0.getFacesKPtr(), 4 * rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_total_ks, g0.getTotalKPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_active, g0.getActivePtr(), rows * cols * sizeof(unsigned char), cudaMemcpyHostToDevice);

 
    int interior = (rows - 2) * (cols - 2); 
    int threadsPerBlock = 256;
    int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
    for (int i = 0; i < sweeps; i++){
        jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols);
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