#ifndef UPDATE_CELL_HPP
#define UPDATE_CELL_HPP

#ifdef __CUDACC__
#define NUMKIT_HD __host__ __device__
#else
#define NUMKIT_HD
#endif 
 
NUMKIT_HD double neighbor_sum(const double* temps, const double* faces_k, int r, int c, int rows, int cols){
    double k_up = faces_k[(r * cols + c) * 4 + 0]; double k_down = faces_k[(r * cols + c) * 4 + 1];
    double k_right = faces_k[(r * cols + c) * 4 + 2]; double k_left = faces_k[(r * cols + c) * 4 + 3];
    return (k_up * temps[(r + 1) * cols + c] + k_down * temps[(r - 1) * cols + c]
      + k_right * temps[r * cols + (c + 1)] + k_left * temps[r * cols + (c - 1)]);
}

NUMKIT_HD double update_cell(const double* temps, const double* faces_k, const double* total_ks, const unsigned char* active, int r, int c, int rows, int cols){
    if (active[r * cols + c]){return neighbor_sum(temps, faces_k, r, c, rows, cols) / total_ks[r * cols + c];}
    else{return temps[r * cols + c];}
}

NUMKIT_HD double neighbor_sum(const double* temps, const double* faces_k, int r, int c, int d, int rows, int cols, int depth){
    double k_up = faces_k[((d * rows * cols) + r * cols + c) * 6 + 0]; double k_down = faces_k[((d * rows * cols) + r * cols + c) * 6 + 1];
    double k_right = faces_k[((d * rows * cols) + r * cols + c) * 6 + 2]; double k_left = faces_k[((d * rows * cols) + r * cols + c) * 6 + 3];
    double k_front = faces_k[((d * rows * cols) + r * cols + c) * 6 + 4]; double k_back = faces_k[((d * rows * cols) + r * cols + c) * 6 + 5];

    return (k_up * temps[(d * rows * cols) + (r + 1) * cols + c] + k_down * temps[(d * rows * cols) + (r - 1) * cols + c]
      + k_right * temps[(d * rows * cols) + r * cols + c + 1] + k_left * temps[(d * rows * cols) + r * cols + c - 1]
      + k_front * temps[((d + 1) * rows * cols) + r * cols + c] + k_back * temps[((d - 1) * rows * cols) + r * cols + c]);
}

NUMKIT_HD double update_cell(const double* temps, const double* faces_k, const double* total_ks, const unsigned char* active, int r, int c, int d, int rows, int cols, int depth){
    if (active[(d * rows * cols) + r * cols + c]){return neighbor_sum(temps, faces_k, r, c, d, rows, cols, depth) / total_ks[(d * rows * cols) + r * cols + c];}
    else{return temps[(d * rows * cols) + r * cols + c];}
}
 
#endif