#include <iostream>
#include <numbers>
#include <cmath>
#include <vector>
#include "grid.hpp"
#include "jacobi.hpp"
#include "jacobi_tiled_kernel.cuh"

double exact(double x, double y){
    return (std::sin(std::numbers::pi * x) * std::sinh(std::numbers::pi * y));
}

void setExactBoundary(Grid& g){
    double h = 1.0/(g.getRows() - 1);
    for (int i = 0; i < g.getRows(); i++){
        g.at(i,0) = exact(i * h, 0);
        g.at(i, g.getCols()-1) = exact(i * h, (g.getCols() - 1) * h);
        g.at(0, i) = exact(0, i * h);
        g.at(g.getRows() - 1, i) = exact((g.getRows() - 1) * h, i * h);
    }
}

double maxError(Grid& computed){
    double maxErr = 0;
    double h = 1.0/(computed.getRows() - 1);
    for (int i = 0; i < computed.getRows(); i++){
        for (int j = 0; j < computed.getCols(); j++){
            double error = std::abs(computed.at(i,j) - exact(i*h,j*h));
            if (error > maxErr){
                maxErr = error;
            }
        }
    }
    return maxErr;
}

int main(){
    std::vector<double> maxErrs;
    for (int k = 20; k <= 160; k *= 2){
        Grid g(k,k);
        Grid g0(k,k);
        setExactBoundary(g);
        for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){
                g0.at(i,j) = g.at(i,j);
            }
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

        dim3 threadsPerBlock(8,8);
        dim3 numBlocks((rows + threadsPerBlock.x - 1) / threadsPerBlock.x, (cols + threadsPerBlock.x - 1) / threadsPerBlock.x);

        for (int i = 0; i < sweeps; i++){
            jacobi_tiled_kernel<<<numBlocks, threadsPerBlock, (threadsPerBlock.x + 2) * (threadsPerBlock.x + 2) * sizeof(double)>>>(d_old, d_new, rows, cols);
            cudaDeviceSynchronize();

            double* temp = d_old;
            d_old = d_new;
            d_new = temp;
        }
        double* g1 = (double*)malloc(rows * cols * sizeof(double));
        cudaMemcpy(g1, d_old, rows * cols * sizeof(double), cudaMemcpyDeviceToHost);
        Grid computed(k,k);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                computed.at(i,j) = g1[i * cols + j];
            }
        }

        double err = maxError(computed);
        maxErrs.push_back(err);
        cudaFree(d_old); cudaFree(d_new); free(g1);
    }
    for (int k = 20; k <= 160; k *= 2){
        std::cout << "Grid Size: " << k << " x " << k << " Error: " << maxErrs[std::log2(k/20)] << "\n";
    }
    return 0;
}