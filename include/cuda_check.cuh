#ifndef CUDA_CHECK_CUH
#define CUDA_CHECK_CUH

#include <cstdio>
#include <cstdlib>

#define CUDA_CHECK(call) do {                                  \
    cudaError_t err = call;                                    \
    if (err != cudaSuccess) {                                  \
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err)); \
        exit(1);                                                \
    }                                                            \
} while (0)

#endif