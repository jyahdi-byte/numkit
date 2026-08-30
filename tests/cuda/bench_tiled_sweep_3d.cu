#include <chrono>
#include <iostream>
#include <vector>
#include <algorithm>
#include "grid.hpp"
#include "jacobi_tiled_kernel.cuh"
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
    double* d_faces_k;
    double* d_total_ks;
    unsigned char* d_active;
    CUDA_CHECK(cudaMalloc((void**)&d_old, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_new, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 6 * n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_total_ks, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&d_active, n * sizeof(unsigned char)));
    CUDA_CHECK(cudaMemcpy(d_old, g.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_new, g.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_faces_k, g.getFacesKPtr(), 6 * n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_total_ks, g.getTotalKPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_active, g.getActivePtr(), n * sizeof(unsigned char), cudaMemcpyHostToDevice));
    std::vector<double> meanTimes;
    std::vector<double> devTimes;
    std::vector<int> blockSizes;

    for (int blockSize = 4; blockSize <= 8; blockSize += 2){
        dim3 threadsPerBlock(blockSize, blockSize, blockSize);
        dim3 numBlocks((rows + threadsPerBlock.x - 1) / threadsPerBlock.x,
                        (cols + threadsPerBlock.y - 1) / threadsPerBlock.y,
                        (depth + threadsPerBlock.z - 1) / threadsPerBlock.z);
        size_t shared = (threadsPerBlock.x + 2) * (threadsPerBlock.y + 2) * (threadsPerBlock.z + 2) * sizeof(double);
        std::vector<double> times;

        for (int i = 0; i < 10; i++){
            jacobi_tiled_kernel<<<numBlocks, threadsPerBlock, shared>>>(d_old, d_new, d_faces_k, d_total_ks, d_active, rows, cols, depth);
            CUDA_CHECK(cudaDeviceSynchronize());

            double* temp = d_old;
            d_old = d_new;
            d_new = temp;
        }

        for (int k = 0; k < 10; k++){

            auto t1 = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < sweeps; i++){
                jacobi_tiled_kernel<<<numBlocks, threadsPerBlock, shared>>>(d_old, d_new, d_faces_k, d_total_ks, d_active, rows, cols, depth);
                CUDA_CHECK(cudaDeviceSynchronize());

                double* temp = d_old;
                d_old = d_new;
                d_new = temp;
            }

            auto t2 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
            std::cout << "Elapsed: " << ms << " ms for " << sweeps << " sweeps, block " << blockSize << "\n";
            times.push_back(ms);
        }
        meanTimes.push_back(mean(times));
        devTimes.push_back(dev(times));
        blockSizes.push_back(blockSize);
        times.clear();
    }
    double min_value = *std::min_element(meanTimes.begin(), meanTimes.end());
    auto min_iterator = std::min_element(meanTimes.begin(), meanTimes.end());
    int min_index = std::distance(meanTimes.begin(), min_iterator);

    std::cout << "\n";
    std::cout << "Mean of Mean Elapsed Times: " << mean(meanTimes) << "\n";
    std::cout << "Mean of Standard Deviation Elapsed Times: " << dev(devTimes) << "\n";
    std::cout << "Fastest Block Size: " << blockSizes[min_index] << " Mean Time: " << min_value << " Standard Deviation: " << devTimes[min_index] << "\n";

    CUDA_CHECK(cudaFree(d_old));
    CUDA_CHECK(cudaFree(d_new));
    CUDA_CHECK(cudaFree(d_faces_k));
    CUDA_CHECK(cudaFree(d_total_ks));
    CUDA_CHECK(cudaFree(d_active));

    return 0;
}