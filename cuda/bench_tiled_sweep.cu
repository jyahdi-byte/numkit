#include <chrono>
#include <iostream>
#include <vector>
#include <algorithm>
#include "grid.hpp"
#include "jacobi_tiled_kernel.cuh"
#include "stats.hpp"

int main(){
    Grid g(100,100);
    for (int j = 0; j < g.getCols(); j++){
        g.at(0,j) = 100;
    }

    int rows = g.getRows();
    int cols = g.getCols();
    int sweeps = 10000;
    double* d_new;
    double* d_old;
    cudaMalloc((void**)&d_new, rows * cols * sizeof(double));
    cudaMalloc((void**)&d_old, rows * cols * sizeof(double));
    cudaMemcpy(d_old, g.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_new, g.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    std::vector<double> meanTimes;
    std::vector<double> devTimes;

    for (int blockSize = 8; blockSize <= 32; blockSize += 8){
        dim3 threadsPerBlock(blockSize, blockSize);
        dim3 numBlocks((rows + threadsPerBlock.x - 1) / threadsPerBlock.x, (cols + threadsPerBlock.y - 1) / threadsPerBlock.y);
        std::vector<double> times;

        for (int i = 0; i < 10; i++){
                jacobi_tiled_kernel<<<numBlocks, threadsPerBlock, (threadsPerBlock.x + 2) * (threadsPerBlock.x + 2) * sizeof(double) >>>(d_old, d_new, rows, cols);
                cudaDeviceSynchronize();

                double* temp = d_old;
                d_old = d_new; 
                d_new = temp;
            }

        for (int k = 0; k < 10; k++){

            // AI added: timer around the real timed loop
            auto t1 = std::chrono::high_resolution_clock::now();

            for (int i = 0; i < sweeps; i++){
                jacobi_tiled_kernel<<<numBlocks, threadsPerBlock, (threadsPerBlock.x + 2) * (threadsPerBlock.x + 2) * sizeof(double) >>>(d_old, d_new, rows, cols);
                cudaDeviceSynchronize();

                double* temp = d_old;
                d_old = d_new;
                d_new = temp;
            }

            // AI added: stop timer, convert to ms
            auto t2 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
            std::cout << "Elapsed: " << ms << " ms for " << sweeps << " sweeps\n";
            times.push_back(ms);
        }
        meanTimes.push_back(mean(times));
        devTimes.push_back(dev(times));
        times.clear();
    }
    double min_value = *std::min_element(meanTimes.begin(), meanTimes.end());
    auto min_iterator = std::min_element(meanTimes.begin(), meanTimes.end());
    int min_index = std::distance(meanTimes.begin(), min_iterator);

    std::cout << "\n";
    std::cout << "Mean of Mean Elapsed Times: " << mean(meanTimes) << "\n";
    std::cout << "Mean of Standard Deviation Elapsed Times: " << dev(devTimes) << "\n";
    std::cout << "Fastest Block Size: " << min_index * 8 + 8 << " Mean Time: " << min_value <<
                 " Standard Deviation: " << devTimes[min_index] << "\n";
    

    return 0;   
}