#include <chrono>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "jacobi_kernel.cuh"
#include "stats.hpp" 
#include "cuda_check.cuh"

int main(){
    Grid g(100,100);
    for (int j = 0; j < g.getCols(); j++){g.at(0,j) = 100;}

    int rows = g.getRows();
    int cols = g.getCols();
    int sweeps = 10000;
    double* d_old;
    double* d_new;
    CellType* d_types;
    double* d_faces_k;
    double* d_total_ks;
    unsigned char* d_active;
    CUDA_CHECK(cudaMalloc((void**)&d_old, rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_new, rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_types, rows * cols * sizeof(CellType)));
    CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 4 * rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_total_ks, rows * cols * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_active, rows * cols * sizeof(unsigned char)));
    CUDA_CHECK(cudaMemcpy(d_old, g.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_new, g.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_types, g.getTypesPtr(), rows * cols * sizeof(CellType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_faces_k, g.getFacesKPtr(), 4 * rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_total_ks, g.getTotalKPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_active, g.getActivePtr(), rows * cols * sizeof(unsigned char), cudaMemcpyHostToDevice));

    int interior = (rows - 2) * (cols - 2);
    int threadsPerBlock = 256;
    int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
    std::vector<double> times;

    for (int i = 0; i < 10; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols);
            CUDA_CHECK(cudaDeviceSynchronize());

            double* temp = d_old;
            d_old = d_new; 
            d_new = temp;
        }

    for (int k = 0; k < 10; k++){

        auto t1 = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < sweeps; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols);
            CUDA_CHECK(cudaDeviceSynchronize());

            double* temp = d_old;
            d_old = d_new;
            d_new = temp;
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        std::cout << "Elapsed: " << ms << " ms for " << sweeps << " sweeps\n";
        times.push_back(ms);
    }
    std::cout << "\n";
    std::cout << "Mean Elapsed Time: " << mean(times) << "\n";
    std::cout << "Elapsed Time Standard Deviation: " << dev(times) << "\n";

    CUDA_CHECK(cudaFree(d_old));
    CUDA_CHECK(cudaFree(d_new));
    CUDA_CHECK(cudaFree(d_types));
    CUDA_CHECK(cudaFree(d_faces_k));
    CUDA_CHECK(cudaFree(d_total_ks));
    CUDA_CHECK(cudaFree(d_active));

    return 0;   
}