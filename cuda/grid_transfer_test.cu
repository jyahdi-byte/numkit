#include <stdio.h>
#include <cassert>
#include "grid.hpp"

int main(){
  int rows = 4;
  int cols = 4;
  Grid g(rows, cols);
  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      g.at(i,j) = i * j;
    }
  }

  double* d_grid;
  cudaMalloc((void**)&d_grid, rows * cols * sizeof(double));
  cudaMemcpy(d_grid, g.getTempsPtr(), rows * cols * sizeof(double), cudaMemcpyHostToDevice);

  double* g1 = (double*)malloc(rows * cols * sizeof(double));
  cudaMemcpy(g1, d_grid, rows * cols * sizeof(double), cudaMemcpyDeviceToHost);

  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      assert(g.at(i,j) == g1[i*cols+j]);
    }
  }

  printf("grid_transfer_test: PASS\n");

  free(g1);
  cudaFree(d_grid);
  return 0;
}
