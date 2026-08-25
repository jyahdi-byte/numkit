#ifndef JACOBI_MT_HPP
#define JACOBI_MT_HPP

#include <thread>
#include <vector>
#include <barrier>
#include <atomic>
#include <functional>
#include <cmath>
#include <algorithm>

#include "grid.hpp" 
#include "update_cell.hpp"

void sweep_rows(const Grid& g, Grid& newg, int r_start, int r_end, double& localMax){
    localMax = 0;
    double change = 0;
    for (int i = r_start; i < r_end; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
            if (g.getType(i,j) == INTERIOR){
                    newg.at(i,j) = update_cell(g.getTempsPtr(), g.getFacesKPtr(), g.getTotalKPtr(), g.getActivePtr(), i, j, g.getRows(), g.getCols());
                    change = std::abs(newg.at(i,j) - g.at(i,j));
                    if (change >= localMax){
                        localMax = change; 
                    }
                }
            }
        } 
}

// persistent version, thread just keeps sweeping its rows until converged flag flips
void worker_job(Grid& g, Grid& newg, int r_start, int r_end, double& localMax, int max_iter, std::barrier<std::function<void()>>& sync_point, std::atomic<bool>& converged){
    for (int k = 0; k < max_iter; k++){
        if (converged.load()) return;
        sweep_rows(g, newg, r_start, r_end, localMax);
        sync_point.arrive_and_wait(); // last thread in triggers on_complete below
    }
}

int jacobi_mt_solve(Grid& g, double tol, int max_iter, int num_threads){
    Grid newg = g;
    int interior = g.getRows() - 2;
    int chunk = interior / num_threads;
    std::vector<double> localMaxes(num_threads, 0.0);
    std::atomic<bool> converged{false};
    int iters_done = 0;

    // swap + tol check, runs once per round instead of every sweep spawning threads
    std::function<void()> on_complete = [&](){
        g = newg;
        double maxChange = *std::max_element(localMaxes.begin(), localMaxes.end());
        iters_done++;
        if (maxChange <= tol) converged.store(true);
    };
    std::barrier<std::function<void()>> sync_point(num_threads, on_complete);

    std::vector<std::thread> workers;
    for (int t = 0; t < num_threads; t++){
        int r_start = 1 + t * chunk;
        int r_end = (t == num_threads - 1) ? g.getRows() - 1 : r_start + chunk;
        workers.push_back(std::thread(worker_job, std::ref(g), std::ref(newg), r_start, r_end, std::ref(localMaxes[t]), max_iter, std::ref(sync_point), std::ref(converged)));
    }
    for (auto& w : workers) w.join();

    return iters_done;
}

#endif