#include <cassert>
#include <cmath>
#include <iostream>

#include "grid1d.hpp"
#include "space_time_log.hpp"
#include "diffusion.hpp"
#include "black_scholes.hpp"
#include "black_scholes_exact.hpp"

int main(){
    int K = 100;
    double sigma = 0.2;
    double rate = 0.05;
    double T = 1;
    double x_min = 2;
    double x_max = 7.5;
    double dx = 0.005;
    double dtau = 0.00025;
    double num_steps = T / dtau;
    double alpha = std::pow(sigma, 2) / 2;
    int points = (int)std::round((x_max - x_min)/dx) + 1;
    
    Grid1D u(points, dx);
    for (int i = 0; i < u.getPoints(); i++){
        double x = x_min + i * dx;
        u.at(i) = std::max(std::exp(x) - K, 0.0);
    }

    black_scholes_solve(u, num_steps, x_max, rate, dtau, sigma, K, alpha, dx);
    int idx = (int)std::round((std::log(100.0) - x_min)/dx);
    double tau_final = num_steps * dtau;
    u.at(idx) = u.at(idx) * std::exp(-rate * tau_final);
    double exact = compute_price(100, K, T, rate, sigma);
    if(std::abs(u.at(idx) - exact) > 0.01 * exact){
        std::cout << "Error of " << std::abs(u.at(idx) - exact) << " is greater than limit of " << 0.01 * exact << std::endl;
        assert(false); 
    }
    std::cout << "PASS\n";
    return 0;
}