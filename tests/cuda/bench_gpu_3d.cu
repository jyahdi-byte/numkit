#include <chrono>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "jacobi_kernel.cuh"
#include "stats.hpp"
#include "cuda_check.cuh"

int main(){
    Grid3D g(40,40,40);
    for (int j = 0; j < g.getCols(); j++){
        for (int k = 0; k < g.getDepth(); k++){
            g.at(0,j,k) = 100;
        }
    }

    int rows = g.getRows();
    int cols = g.getCols();
    int depth = g.getDepth();
    int n = rows * cols * depth;
    int sweeps = 5000;
    double* d_old;
    double* d_new;
    CellType* d_types;
    double* d_faces_k;
    double* d_total_ks;
    unsigned char* d_active;
    CUDA_CHECK(cudaMalloc((void**)&d_old, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_new, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_types, n * sizeof(CellType)));
    CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 6 * n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_total_ks, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_active, n * sizeof(unsigned char)));
    CUDA_CHECK(cudaMemcpy(d_old, g.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_new, g.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_types, g.getTypesPtr(), n * sizeof(CellType), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_faces_k, g.getFacesKPtr(), 6 * n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_total_ks, g.getTotalKPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_active, g.getActivePtr(), n * sizeof(unsigned char), cudaMemcpyHostToDevice));

    int interior = (rows - 2) * (cols - 2) * (depth - 2);
    int threadsPerBlock = 256;
    int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
    std::vector<double> times;

    for (int i = 0; i < 10; i++){
        jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth);
        CUDA_CHECK(cudaDeviceSynchronize());

        double* temp = d_old;
        d_old = d_new;
        d_new = temp;
    }

    for (int k = 0; k < 10; k++){

        auto t1 = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < sweeps; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth);
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