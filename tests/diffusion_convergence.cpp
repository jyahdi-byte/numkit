#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

#include "grid1d.hpp"
#include "space_time_log.hpp"
#include "diffusion.hpp"
#include "diffusion_exact.hpp"

int main(){

    std::vector<double> dx_vector;
    double dx0 = 0.04;
    for (int i = 0; i < 5; i++){
        dx_vector.push_back(dx0 / (std::pow(2, i)));
    }
    double A0 = 1;
    double alpha = 1;
    double num_steps = 125; 
    double dt = 0.0004;
    std::vector<double> expected_error_vector = {3 * 3.05e-04, 3 * 7.61e-05, 3 * 1.90e-05, 3 * 4.75e-06, 3 * 1.86e-06};

    for (size_t i = 0; i < dx_vector.size(); i++){
        std::cout << "Test for dx = " << dx_vector[i] << "\n";
        double sigma = 6*dx_vector[0];
        Grid1D u(100 * std::pow(2,i), dx_vector[i]);
        double x0 = u.getLength()/2;
        for (int j = 0; j < u.getPoints(); j++){
            u.at(j) = compute_u_exact(j * dx_vector[i], 0, x0, sigma, A0, alpha);
        }

        num_steps *= 4;
        dt /= 4;
        double t_actual = num_steps * dt;
        double r = compute_r(alpha, dt, dx_vector[i]);
        SpaceTimeLog log(num_steps + 1, u.getPoints());
        diffusion_solve(u, num_steps, r, log);
        double maxErr = 0;  
        for (int j = 0; j < u.getPoints(); j++){
            double diff = std::abs(compute_u_exact(j * dx_vector[i], t_actual, x0, sigma, A0, alpha) - u.at(j));
            if (diff > maxErr){maxErr = diff;}
        }

        if (maxErr > expected_error_vector[i]){
            std::cout << "Max Error: " << maxErr << ", is greater than limit of " << expected_error_vector[i] << std::endl;
            assert(false);
        }

        std::cout << "Max Error: " << maxErr << "\n";
        std::cout << "PASS\n\n";
}

    return 0;
}