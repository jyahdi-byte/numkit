#include <iostream>
#include <cassert>
#include <cmath>
#include <numbers>
#include "grid.hpp"
#include "jacobi.hpp"
#include "gauss_seidel_rb.hpp"
#include "sor_rb.hpp"
#include "jacobi_kernel.cuh"
#include "gauss_seidel_rb_kernel.cuh"
#include "sor_rb_kernel.cuh"
#include "jacobi_tiled_kernel.cuh"
#include "cuda_check.cuh"

Grid3D make_grid(){
    Grid3D g(10,10,10);
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            for (int k = 0; k < g.getDepth(); k++){
                if (i == 0 || j == 0 || k == 0 || i == g.getRows()-1 || j == g.getCols()-1 || k == g.getDepth()-1){
                    g.at(i,j,k) = i;
                }
            }
        }
    }
    g.maskRect(4,4,4,6,6,6,HOLE);
    return g;
}

int main(){

    // ============================================================
    // 1. Jacobi (untiled) vs CPU jacobi_solve
    // ============================================================
    {
        Grid3D g = make_grid();
        Grid3D g0 = g;

        int rows = g.getRows();
        int cols = g.getCols();
        int depth = g.getDepth();
        int n = rows * cols * depth;

        int sweeps = jacobi_solve(g, 1e-10, 10000);
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
        CUDA_CHECK(cudaMemcpy(d_old, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_new, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_types, g0.getTypesPtr(), n * sizeof(CellType), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_faces_k, g0.getFacesKPtr(), 6 * n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_total_ks, g0.getTotalKPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_active, g0.getActivePtr(), n * sizeof(unsigned char), cudaMemcpyHostToDevice));

        int interior = (rows - 2) * (cols - 2) * (depth - 2);
        int threadsPerBlock = 256;
        int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
        for (int i = 0; i < sweeps; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(d_old, d_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth);
            CUDA_CHECK(cudaDeviceSynchronize());

            double* temp = d_old;
            d_old = d_new;
            d_new = temp;
        }
        double* g1 = (double*)malloc(n * sizeof(double));
        CUDA_CHECK(cudaMemcpy(g1, d_old, n * sizeof(double), cudaMemcpyDeviceToHost));

        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                for (int k = 0; k < depth; k++){
                    assert(std::abs(g.at(i,j,k) - g1[(k * rows * cols) + i * cols + j]) < 1e-9);
                }
            }
        }
        std::cout << "PASS: jacobi_kernel (3D), " << sweeps << " sweeps\n";

        free(g1);
        CUDA_CHECK(cudaFree(d_old));
        CUDA_CHECK(cudaFree(d_new));
        CUDA_CHECK(cudaFree(d_types));
        CUDA_CHECK(cudaFree(d_faces_k));
        CUDA_CHECK(cudaFree(d_total_ks));
        CUDA_CHECK(cudaFree(d_active));
    }

    // ============================================================
    // 2. Gauss-Seidel red-black vs CPU gauss_seidel_rb_solve
    // ============================================================
    {
        Grid3D g = make_grid();
        Grid3D g0 = g;

        int rows = g.getRows();
        int cols = g.getCols();
        int depth = g.getDepth();
        int n = rows * cols * depth;

        int sweeps = gauss_seidel_rb_solve(g, 1e-10, 10000);
        double* d;
        CellType* d_types;
        double* d_faces_k;
        double* d_total_ks;
        unsigned char* d_active;

        CUDA_CHECK(cudaMalloc((void**)&d, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_types, n * sizeof(CellType)));
        CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 6 * n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_total_ks, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_active, n * sizeof(unsigned char)));
        CUDA_CHECK(cudaMemcpy(d, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_types, g0.getTypesPtr(), n * sizeof(CellType), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_faces_k, g0.getFacesKPtr(), 6 * n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_total_ks, g0.getTotalKPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_active, g0.getActivePtr(), n * sizeof(unsigned char), cudaMemcpyHostToDevice));

        int interior = (rows - 2) * (cols - 2) * (depth - 2);
        int threadsPerBlock = 256;
        int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
        for (int i = 0; i < sweeps; i++){
            gs_rb_kernel_r<<<numBlocks, threadsPerBlock>>>(d, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth);
            CUDA_CHECK(cudaDeviceSynchronize());
            gs_rb_kernel_b<<<numBlocks, threadsPerBlock>>>(d, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth);
            CUDA_CHECK(cudaDeviceSynchronize());
        }
        double* g1 = (double*)malloc(n * sizeof(double));
        CUDA_CHECK(cudaMemcpy(g1, d, n * sizeof(double), cudaMemcpyDeviceToHost));

        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                for (int k = 0; k < depth; k++){
                    assert(std::abs(g.at(i,j,k) - g1[(k * rows * cols) + i * cols + j]) < 1e-9);
                }
            }
        }
        std::cout << "PASS: gs_rb_kernel (3D), " << sweeps << " sweeps\n";

        free(g1);
        CUDA_CHECK(cudaFree(d));
        CUDA_CHECK(cudaFree(d_types));
        CUDA_CHECK(cudaFree(d_faces_k));
        CUDA_CHECK(cudaFree(d_total_ks));
        CUDA_CHECK(cudaFree(d_active));
    }

    // ============================================================
    // 3. SOR red-black vs CPU sor_rb_solve
    // ============================================================
    {
        Grid3D g = make_grid();
        Grid3D g0 = g;

        int rows = g.getRows();
        int cols = g.getCols();
        int depth = g.getDepth();
        int n = rows * cols * depth;
        double h = 1.0 / (rows - 1);
        double omega = 2.0 / (1.0 + std::sin(std::numbers::pi * h));

        int sweeps = sor_rb_solve(g, 1e-10, 10000, omega);
        double* d;
        CellType* d_types;
        double* d_faces_k;
        double* d_total_ks;
        unsigned char* d_active;

        CUDA_CHECK(cudaMalloc((void**)&d, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_types, n * sizeof(CellType)));
        CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 6 * n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_total_ks, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_active, n * sizeof(unsigned char)));
        CUDA_CHECK(cudaMemcpy(d, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_types, g0.getTypesPtr(), n * sizeof(CellType), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_faces_k, g0.getFacesKPtr(), 6 * n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_total_ks, g0.getTotalKPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_active, g0.getActivePtr(), n * sizeof(unsigned char), cudaMemcpyHostToDevice));

        int interior = (rows - 2) * (cols - 2) * (depth - 2);
        int threadsPerBlock = 256;
        int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
        for (int i = 0; i < sweeps; i++){
            sor_rb_kernel_r<<<numBlocks, threadsPerBlock>>>(d, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth, omega);
            CUDA_CHECK(cudaDeviceSynchronize());
            sor_rb_kernel_b<<<numBlocks, threadsPerBlock>>>(d, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth, omega);
            CUDA_CHECK(cudaDeviceSynchronize());
        }
        double* g1 = (double*)malloc(n * sizeof(double));
        CUDA_CHECK(cudaMemcpy(g1, d, n * sizeof(double), cudaMemcpyDeviceToHost));

        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                for (int k = 0; k < depth; k++){
                    assert(std::abs(g.at(i,j,k) - g1[(k * rows * cols) + i * cols + j]) < 1e-9);
                }
            }
        }
        std::cout << "PASS: sor_rb_kernel (3D), " << sweeps << " sweeps\n";

        free(g1);
        CUDA_CHECK(cudaFree(d));
        CUDA_CHECK(cudaFree(d_types));
        CUDA_CHECK(cudaFree(d_faces_k));
        CUDA_CHECK(cudaFree(d_total_ks));
        CUDA_CHECK(cudaFree(d_active));
    }

    // ============================================================
    // 4. Tiled Jacobi vs untiled Jacobi kernel (GPU vs GPU, per step 5)
    // ============================================================
    {
        Grid3D g = make_grid();
        Grid3D g0 = g;

        int rows = g.getRows();
        int cols = g.getCols();
        int depth = g.getDepth();
        int n = rows * cols * depth;
        int sweeps = 5000;

        CellType* d_types;
        double* d_faces_k;
        double* d_total_ks;
        unsigned char* d_active;
        CUDA_CHECK(cudaMalloc((void**)&d_types, n * sizeof(CellType)));
        CUDA_CHECK(cudaMalloc((void**)&d_faces_k, 6 * n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_total_ks, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_active, n * sizeof(unsigned char)));
        CUDA_CHECK(cudaMemcpy(d_types, g0.getTypesPtr(), n * sizeof(CellType), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_faces_k, g0.getFacesKPtr(), 6 * n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_total_ks, g0.getTotalKPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_active, g0.getActivePtr(), n * sizeof(unsigned char), cudaMemcpyHostToDevice));

        // untiled reference run
        double* r_old;
        double* r_new;
        CUDA_CHECK(cudaMalloc((void**)&r_old, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&r_new, n * sizeof(double)));
        CUDA_CHECK(cudaMemcpy(r_old, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(r_new, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));

        int interior = (rows - 2) * (cols - 2) * (depth - 2);
        int threadsPerBlock = 256;
        int numBlocks = (interior + threadsPerBlock - 1) / threadsPerBlock;
        for (int i = 0; i < sweeps; i++){
            jacobi_kernel<<<numBlocks, threadsPerBlock>>>(r_old, r_new, d_types, d_faces_k, d_total_ks, d_active, rows, cols, depth);
            CUDA_CHECK(cudaDeviceSynchronize());

            double* temp = r_old;
            r_old = r_new;
            r_new = temp;
        }
        double* g_ref = (double*)malloc(n * sizeof(double));
        CUDA_CHECK(cudaMemcpy(g_ref, r_old, n * sizeof(double), cudaMemcpyDeviceToHost));

        // tiled run
        double* d_old;
        double* d_new;
        CUDA_CHECK(cudaMalloc((void**)&d_old, n * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&d_new, n * sizeof(double)));
        CUDA_CHECK(cudaMemcpy(d_old, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_new, g0.getTempsPtr(), n * sizeof(double), cudaMemcpyHostToDevice));

        dim3 threadsPerBlock3(8,8,8);
        dim3 numBlocks3((rows + threadsPerBlock3.x - 1) / threadsPerBlock3.x,
                         (cols + threadsPerBlock3.y - 1) / threadsPerBlock3.y,
                         (depth + threadsPerBlock3.z - 1) / threadsPerBlock3.z);

        for (int i = 0; i < sweeps; i++){
            jacobi_tiled_kernel<<<numBlocks3, threadsPerBlock3, (threadsPerBlock3.x + 2) * (threadsPerBlock3.y + 2) * (threadsPerBlock3.z + 2) * sizeof(double)>>>(d_old, d_new, d_faces_k, d_total_ks, d_active, rows, cols, depth);
            CUDA_CHECK(cudaDeviceSynchronize());

            double* temp = d_old;
            d_old = d_new;
            d_new = temp;
        }
        double* g1 = (double*)malloc(n * sizeof(double));
        CUDA_CHECK(cudaMemcpy(g1, d_old, n * sizeof(double), cudaMemcpyDeviceToHost));

        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                for (int k = 0; k < depth; k++){
                    int idx = (k * rows * cols) + i * cols + j;
                    assert(std::abs(g_ref[idx] - g1[idx]) < 1e-9);
                }
            }
        }
        std::cout << "PASS: jacobi_tiled_kernel matches untiled (3D), " << sweeps << " sweeps\n";

        free(g_ref);
        free(g1);
        CUDA_CHECK(cudaFree(d_types));
        CUDA_CHECK(cudaFree(d_faces_k));
        CUDA_CHECK(cudaFree(d_total_ks));
        CUDA_CHECK(cudaFree(d_active));
        CUDA_CHECK(cudaFree(r_old));
        CUDA_CHECK(cudaFree(r_new));
        CUDA_CHECK(cudaFree(d_old));
        CUDA_CHECK(cudaFree(d_new));
    }

    std::cout << "\nALL 3D CUDA KERNELS PASS\n";
    return 0;
}