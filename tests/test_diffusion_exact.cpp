#include <cassert>
#include <iostream>
#include <cmath>

#include "grid1d.hpp"
#include "space_time_log.hpp"
#include "diffusion.hpp"
#include "ppm.hpp"

double compute_sigmaT(double sigma, double alpha, double t){
    return std::sqrt(std::pow(sigma,2) + 2 * alpha * t);
}

double compute_At(double A0, double sigma, double alpha, double t){
    double sigmaT = compute_sigmaT(sigma, alpha, t);
    return A0 * sigma / sigmaT; 
}

double compute_u_exact(double x, double t, double x0, double sigma, double A0, double alpha){
    double sigmaT = compute_sigmaT(sigma, alpha, t);
    double aT = compute_At(A0, sigma, alpha, t);
    return aT * std::exp(-std::pow((x-x0 ),2) / (2 * std::pow(sigmaT,2)));
}

int main(){

    double dx = 0.02;
    double A0 = 1;
    double alpha = 1;
    double sigma = 6*dx;
    Grid1D u(200, dx);
    double x0 = u.getLength()/2;
    for (int i = 0; i < u.getPoints(); i++){
        u.at(i) = compute_u_exact(i * dx, 0, x0, sigma, A0, alpha);
    }

    double num_steps = 500; 
    double dt = 0.0001;
    double t_actual = num_steps * dt;
    double r = compute_r(alpha, dt, dx);
    SpaceTimeLog log(num_steps + 1, u.getPoints());
    diffusion_solve(u, num_steps, r, log);
    double maxErr = 0;  
    for (int i = 0; i < u.getPoints(); i++){
        double diff = std::abs(compute_u_exact(i * dx, t_actual, x0, sigma, A0, alpha) - u.at(i));
        if (diff > maxErr){maxErr = diff;}
    }

    if (maxErr > 2e-4){
        std::cout << "Max Error: " << maxErr << ", is greater than limit of 2e-4, which is ~3x the measured truncation error at this resolution" << std::endl;
        assert(false);
    }

    std::cout << "Max Error: " << maxErr << "\n";
    std::cout << "PASS\n";

    write_ppm(log, "test_diffusion_exact.ppm");

    return 0;
}