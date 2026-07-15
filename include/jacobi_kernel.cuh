#ifndef JACOBI_KERNEL_CUH
#define JACOBI_KERNEL_CUH

__global__ void jacobi_kernel(double* d_old, double* d_new, int rows, int cols) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (rows - 2) * (cols - 2)) return;

    int row = index / (cols - 2);
    int col = index % (cols - 2);
    int i = row + 1;
    int j = col + 1;

    d_new[i * cols + j] = (d_old[i * cols + j - 1] + d_old[i * cols + j + 1] + d_old[(i - 1) * cols + j] + d_old[(i + 1) * cols + j]) / 4;
}

#endif