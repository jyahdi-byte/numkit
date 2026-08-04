#ifndef DIFFUSION_EXACT_HPP
#define DIFFUSION_EXACT_HPP

#include <cmath>

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

#endif
