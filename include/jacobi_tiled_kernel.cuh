#ifndef JACOBI_TILED_KERNEL_CUH
#define JACOBI_TILED_KERNEL_CUH

__global__ void jacobi_tiled_kernel(double* d_old, double* d_new, int rows, int cols) {
    int tileWidth = blockDim.x + 2;
    extern __shared__ double tile[];
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (i < rows-1 && j < cols-1){
        tile[(threadIdx.x + 1) * tileWidth + (threadIdx.y + 1)] = d_old[i * cols + j];
        if (threadIdx.y == 0){ tile[(threadIdx.x + 1) * tileWidth + threadIdx.y] = d_old[i * cols + (j-1)]; }
        if (threadIdx.y == blockDim.x - 1 || j == cols - 2){ tile[(threadIdx.x + 1) * tileWidth + threadIdx.y + 2] = d_old[i * cols + (j+1)]; }
        if (threadIdx.x == 0){ tile[(threadIdx.x) * tileWidth + (threadIdx.y + 1)] = d_old[(i - 1) * cols + j]; }
        if (threadIdx.x == blockDim.x - 1 || i == rows - 2){ tile[(threadIdx.x + 2) * tileWidth + (threadIdx.y + 1)] = d_old[(i + 1) * cols + j]; }
    }   
    __syncthreads(); 

    if ( i < rows-1 && j < cols-1){
        d_new[i * cols + j] = (tile[(threadIdx.x+2) * tileWidth + threadIdx.y + 1] + tile[(threadIdx.x) * tileWidth + threadIdx.y + 1] +
        tile[(threadIdx.x+1) * tileWidth + threadIdx.y+2] + tile[(threadIdx.x+1) * tileWidth + threadIdx.y]) / 4;
    }
}

#endif 