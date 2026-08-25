#ifndef UPDATE_CELL_HPP
#define UPDATE_CELL_HPP

#ifdef __CUDACC__
#define NUMKIT_HD __host__ __device__
#else
#define NUMKIT_HD
#endif

NUMKIT_HD double update_cell(const double* temps, const double* faces_k, const double* total_ks, const unsigned char* active, int r, int c, int rows, int cols){
    if (active[r * cols + c]){
        double k_up = faces_k[(r * cols + c) * 4 + 0];
        double k_down = faces_k[(r * cols + c) * 4 + 1];
        double k_right = faces_k[(r * cols + c) * 4 + 2];
        double k_left = faces_k[(r * cols + c) * 4 + 3];
        return (k_up * temps[(r + 1) * cols + c] + k_down * temps[(r - 1) * cols + c] + k_right * temps[r * cols + (c + 1)] + k_left * temps[r * cols + (c - 1)]) / total_ks[r * cols + c];
    }
    else{return temps[r * cols + c];}
}

#endif