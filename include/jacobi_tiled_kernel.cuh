#ifndef JACOBI_TILED_KERNEL_CUH
#define JACOBI_TILED_KERNEL_CUH

__global__ void jacobi_tiled_kernel(double* d_old, double* d_new, int rows, int cols) {
    __shared__ double tile[18 * 18];
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (i < rows-1 && j < cols-1){
        tile[(threadIdx.x + 1) * 18 + (threadIdx.y + 1)] = d_old[i * cols + j];
        if (threadIdx.y == 0){ tile[(threadIdx.x + 1) * 18 + threadIdx.y] = d_old[i * cols + (j-1)]; }
        if (threadIdx.y == 15){ tile[(threadIdx.x + 1) * 18 + threadIdx.y + 2] = d_old[i * cols + (j+1)]; }
        if (threadIdx.x == 0){ tile[(threadIdx.x) * 18 + (threadIdx.y + 1)] = d_old[(i - 1) * cols + j]; }
        if (threadIdx.x == 15){ tile[(threadIdx.x + 2) * 18 + (threadIdx.y + 1)] = d_old[(i + 1) * cols + j]; }
    }   
    __syncthreads(); 

    if ( i < rows-1 && j < cols-1){
        d_new[i * cols + j] = (tile[(threadIdx.x+2) * 18 + threadIdx.y + 1] + tile[(threadIdx.x) * 18 + threadIdx.y + 1] +
        tile[(threadIdx.x+1) * 18 + threadIdx.y+2] + tile[(threadIdx.x+1) * 18 + threadIdx.y]) / 4;
    }
}

#endif