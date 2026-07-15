#include <chrono>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "jacobi_kernel.cuh"
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

    int interior = (rows - 2) * (cols - 2);
    int threadsPerBlock = 256;
    int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
    std::vector<double> times;

    for (int i = 0; i < 10; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, rows, cols);
            cudaDeviceSynchronize();

            double* temp = d_old;
            d_old = d_new;
            d_new = temp;
        }

    for (int k = 0; k < 10; k++){

        // AI added: timer around the real timed loop
        auto t1 = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < sweeps; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, rows, cols);
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
    std::cout << "\n";
    std::cout << "Mean Elapsed Time: " << mean(times) << "\n";
    std::cout << "Elapsed Time Standard Deviation: " << dev(times) << "\n";

    return 0;   
}