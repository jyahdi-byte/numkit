#ifndef GAUSS_SEIDEL_RB_MT_HPP
#define GAUSS_SEIDEL_RB_MT_HPP

#include <thread>
#include <vector>
#include <barrier>
#include <atomic>
#include <functional>
#include <cmath>
#include <algorithm>

#include "grid.hpp"
#include "update_cell.hpp"

void gs_rb_sweep_rows_r(Grid& g, int r_start, int r_end, double& localMax){
    localMax = 0;
    double change = 0;
    for (int i = r_start; i < r_end; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double oldPoint = g.at(i,j);
                if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 0){
                        g.at(i,j) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                        change = std::abs(g.at(i,j) - oldPoint);
                        if (change >= localMax){
                            localMax = change; 
                        }
                    }
            }
        } 
}

void gs_rb_sweep_rows_b(Grid& g, int r_start, int r_end, double& localMax){
    localMax = 0;
    double change = 0;
    for (int i = r_start; i < r_end; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double oldPoint = g.at(i,j);
                if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 1){
                        g.at(i,j) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                        change = std::abs(g.at(i,j) - oldPoint);
                        if (change >= localMax){
                            localMax = change;  
                        }
                    }
            }
        } 
}

void gs_rb_worker_job(Grid& g, int r_start, int r_end, double& localMax, int max_iter, std::barrier<std::function<void()>>& sync_point, std::atomic<bool>& converged){
    for (int k = 0; k < max_iter; k++){
        if (converged.load()) return;
        gs_rb_sweep_rows_r(g, r_start, r_end, localMax);
        sync_point.arrive_and_wait(); 
        gs_rb_sweep_rows_b(g, r_start, r_end, localMax);
        sync_point.arrive_and_wait(); 
    }
}

int gauss_seidel_rb_mt_solve(Grid& g, double tol, int max_iter, int num_threads){
    int interior = g.getRows() - 2;
    int chunk = interior / num_threads;
    std::vector<double> localMaxes(num_threads, 0.0);
    std::atomic<bool> converged{false};
    int iters_done = 0;

    int phase = 0;
    double maxChange_r = 0;
    double maxChange_b = 0;
    std::function<void()> on_complete = [&](){
        phase++;
        if (phase % 2 == 1){maxChange_r = *std::max_element(localMaxes.begin(), localMaxes.end());}
        else if (phase % 2 == 0){
            maxChange_b = *std::max_element(localMaxes.begin(), localMaxes.end());
            iters_done++;
            double maxChange = std::max(maxChange_r, maxChange_b);
            if (maxChange <= tol){ converged.store(true);}
        }
    };
    std::barrier<std::function<void()>> sync_point(num_threads, on_complete);

    std::vector<std::thread> workers;
    for (int t = 0; t < num_threads; t++){
        int r_start = 1 + t * chunk;
        int r_end = (t == num_threads - 1) ? g.getRows() - 1 : r_start + chunk;
        workers.push_back(std::thread(gs_rb_worker_job, std::ref(g), r_start, r_end, std::ref(localMaxes[t]), max_iter, std::ref(sync_point), std::ref(converged)));
    }
    for (auto& w : workers) w.join();

    return iters_done;
}

void gs_rb_sweep_rows_r(Grid3D& g, int r_start, int r_end, double& localMax){
    localMax = 0;
    double change = 0;
    for (int k = r_start; k < r_end; k++){
        for (int i = 1; i < g.getRows() - 1; i++){
                for (int j = 1; j < g.getCols() - 1; j++){
                    double oldPoint = g.at(i,j,k);
                    if (g.getType(i,j,k) == INTERIOR && (i + j + k) % 2 == 0){
                            g.at(i,j,k) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, k, g.getRows(), g.getCols(), g.getDepth());
                            change = std::abs(g.at(i,j,k) - oldPoint);
                            if (change >= localMax){
                                localMax = change; 
                            }
                        }
                }
            } 
        }
}

void gs_rb_sweep_rows_b(Grid3D& g, int r_start, int r_end, double& localMax){
    localMax = 0;
    double change = 0;
    for (int k = r_start; k < r_end; k++){
        for (int i = 1; i < g.getRows() - 1; i++){
                for (int j = 1; j < g.getCols() - 1; j++){
                    double oldPoint = g.at(i,j,k);
                    if (g.getType(i,j,k) == INTERIOR && (i + j + k) % 2 == 1){
                            g.at(i,j,k) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, k, g.getRows(), g.getCols(), g.getDepth());
                            change = std::abs(g.at(i,j,k) - oldPoint);
                            if (change >= localMax){
                                localMax = change; 
                            }
                        }
                }
            } 
        }
}

void gs_rb_worker_job_3d(Grid3D& g, int r_start, int r_end, double& localMax, int max_iter, std::barrier<std::function<void()>>& sync_point, std::atomic<bool>& converged){
    for (int k = 0; k < max_iter; k++){
        if (converged.load()) return;
        gs_rb_sweep_rows_r(g, r_start, r_end, localMax);
        sync_point.arrive_and_wait(); 
        gs_rb_sweep_rows_b(g, r_start, r_end, localMax);
        sync_point.arrive_and_wait(); 
    }
}

int gauss_seidel_rb_mt_solve(Grid3D& g, double tol, int max_iter, int num_threads){
    int interior = g.getDepth() - 2;
    int chunk = interior / num_threads;
    std::vector<double> localMaxes(num_threads, 0.0);
    std::atomic<bool> converged{false};
    int iters_done = 0;

    int phase = 0;
    double maxChange_r = 0;
    double maxChange_b = 0;
    std::function<void()> on_complete = [&](){
        phase++;
        if (phase % 2 == 1){maxChange_r = *std::max_element(localMaxes.begin(), localMaxes.end());}
        else if (phase % 2 == 0){
            maxChange_b = *std::max_element(localMaxes.begin(), localMaxes.end());
            iters_done++;
            double maxChange = std::max(maxChange_r, maxChange_b);
            if (maxChange <= tol){ converged.store(true);}
        }
    };
    std::barrier<std::function<void()>> sync_point(num_threads, on_complete);

    std::vector<std::thread> workers;
    for (int t = 0; t < num_threads; t++){
        int r_start = 1 + t * chunk;
        int r_end = (t == num_threads - 1) ? g.getDepth() - 1 : r_start + chunk;
        workers.push_back(std::thread(gs_rb_worker_job_3d, std::ref(g), r_start, r_end, std::ref(localMaxes[t]), max_iter, std::ref(sync_point), std::ref(converged)));
    }
    for (auto& w : workers) w.join();

    return iters_done;
}

#endif