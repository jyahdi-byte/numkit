#ifndef SOR_RB_MT_HPP
#define SOR_RB_MT_HPP

#include <thread>
#include <vector>
#include <barrier>
#include <atomic>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numbers>

#include "grid.hpp"

void sor_rb_sweep_rows_r(Grid& g, int r_start, int r_end, double& localMax, double omega){
    localMax = 0;
    double change = 0;
    for (int i = r_start; i < r_end; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double oldPoint = g.at(i,j);
                if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 0){
                        std::vector<double> components = {g.at(i+1,j), g.at(i-1,j), g.at(i,j+1), g.at(i,j-1)}; 
                        if (g.getType(i+1,j) == HOLE){components[0] = g.at(i,j);}
                        if (g.getType(i-1,j) == HOLE){components[1] = g.at(i,j);}
                        if (g.getType(i,j+1) == HOLE){components[2] = g.at(i,j);}
                        if (g.getType(i,j-1) == HOLE){components[3] = g.at(i,j);}
                        double avg = (components[0] + components[1] + components[2] + components[3])/4;
                        double delta = omega * (avg - oldPoint);
                        g.at(i,j) = oldPoint + delta;
                        change = std::abs(delta);
                        if (change >= localMax){
                            localMax = change; 
                        }
                    }
            }
        } 
}

void sor_rb_sweep_rows_b(Grid& g, int r_start, int r_end, double& localMax, double omega){
    localMax = 0;
    double change = 0;
    for (int i = r_start; i < r_end; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double oldPoint = g.at(i,j);
                if (g.getType(i,j) == INTERIOR && (i + j) % 2 == 1){
                        std::vector<double> components = {g.at(i+1,j), g.at(i-1,j), g.at(i,j+1), g.at(i,j-1)}; 
                        if (g.getType(i+1,j) == HOLE){components[0] = g.at(i,j);}
                        if (g.getType(i-1,j) == HOLE){components[1] = g.at(i,j);}
                        if (g.getType(i,j+1) == HOLE){components[2] = g.at(i,j);}
                        if (g.getType(i,j-1) == HOLE){components[3] = g.at(i,j);}
                        double avg = (components[0] + components[1] + components[2] + components[3])/4;
                        double delta = omega * (avg - oldPoint);
                        g.at(i,j) = oldPoint + delta;
                        change = std::abs(delta);
                        if (change >= localMax){
                            localMax = change;  
                        }
                    }
            }
        } 
}

// persistent version, thread just keeps sweeping its rows until converged flag flips
void sor_rb_worker_job(Grid& g, int r_start, int r_end, double& localMax, int max_iter, std::barrier<std::function<void()>>& sync_point, std::atomic<bool>& converged, double omega){
    for (int k = 0; k < max_iter; k++){
        if (converged.load()) return;
        sor_rb_sweep_rows_r(g, r_start, r_end, localMax, omega);
        sync_point.arrive_and_wait(); 
        sor_rb_sweep_rows_b(g, r_start, r_end, localMax, omega);
        sync_point.arrive_and_wait(); 
    }
}

int sor_rb_mt_solve(Grid& g, double tol, int max_iter, int num_threads, double omega){
    int interior = g.getRows() - 2;
    int chunk = interior / num_threads;
    std::vector<double> localMaxes(num_threads, 0.0);
    std::atomic<bool> converged{false};
    int iters_done = 0;

    // swap + tol check, runs once per round instead of every sweep spawning threads
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
        workers.push_back(std::thread(sor_rb_worker_job, std::ref(g), r_start, r_end, std::ref(localMaxes[t]), max_iter, std::ref(sync_point), std::ref(converged), omega));
    }
    for (auto& w : workers) w.join();

    return iters_done;
}

int sor_rb_mt_solve(Grid& g, double tol, int max_iter, int num_threads){
    double h = 1.0 / (g.getRows() - 1);
    double omega = 2.0 / (1.0 + std::sin(std::numbers::pi * h));
    return sor_rb_mt_solve(g, tol, max_iter, num_threads, omega);
}

#endif