#ifndef JACOBI_MT_HPP
#define JACOBI_MT_HPP

#include <thread>
#include <vector>
#include <cmath>
#include "grid.hpp"

double max_element(const std::vector<double> nums){
    double max = nums[0];
    for (size_t i = 0; i < nums.size(); i++){
        if (nums[i] > max){
            max = nums[i];
        }
    }
    return max;
}

void sweep_rows(const Grid& g, Grid& newg, int r_start, int r_end, double& localMax){
    localMax = 0;
    double change = 0;
    for (int i = r_start; i < r_end; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                newg.at(i,j) = (g.at(i+1,j) + g.at(i-1,j) + g.at(i,j+1) + g.at(i,j-1))/4;
                change = std::abs(newg.at(i,j) - g.at(i,j));
                if (change >= localMax){
                    localMax = change;
                }
            }
        }
}

int jacobi_mt_solve(Grid& g, double tol, int max_iter, int num_threads){
    Grid newg = g;
    int interior = g.getRows() - 2;
    int chunk = interior / num_threads;
    
    for (int k = 0; k < max_iter; k++){

        std::vector<double> localMaxes(num_threads, 0.0);
        std::vector<std::thread> workers;
        for (int t = 0; t < num_threads; t++){
            int r_start = 1 + t * chunk;
            int r_end;
            if (t == num_threads - 1){
                r_end = g.getRows() - 1;
            }
            else{
                r_end = r_start + chunk;
            }
            workers.push_back(std::thread(sweep_rows, std::cref(g), std::ref(newg), r_start, r_end, std::ref(localMaxes[t])));
        }
        for(auto& w : workers){
            w.join();
        }
        double maxChange = max_element(localMaxes);
        
        g = newg;
        if (maxChange <= tol){
            return k+1;
        }
    }
    return max_iter;
}

#endif