#ifndef JACOBI_TILED_KERNEL_CUH
#define JACOBI_TILED_KERNEL_CUH

__global__ void jacobi_tiled_kernel(double* d_old, double* d_new, double* d_faces_k, double* d_total_ks, unsigned char* d_active, int rows, int cols) {
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
        int idx = i * cols + j;
        if (d_active[idx]){
            double k_up    = d_faces_k[idx * 4 + 0];
            double k_down  = d_faces_k[idx * 4 + 1];
            double k_right = d_faces_k[idx * 4 + 2];
            double k_left  = d_faces_k[idx * 4 + 3];
            double t_up    = tile[(threadIdx.x+2) * tileWidth + threadIdx.y + 1];
            double t_down  = tile[(threadIdx.x)   * tileWidth + threadIdx.y + 1];
            double t_right = tile[(threadIdx.x+1) * tileWidth + threadIdx.y + 2];
            double t_left  = tile[(threadIdx.x+1) * tileWidth + threadIdx.y];
            d_new[idx] = (k_up*t_up + k_down*t_down + k_right*t_right + k_left*t_left) / d_total_ks[idx];
        } else {
            d_new[idx] = d_old[idx];
        }
    }
}

#endif