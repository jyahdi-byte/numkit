#ifndef BLACK_SCHOLES_EXACT_HPP
#define BLACK_SCHOLES_EXACT_HPP

#include <cmath>

double normal_cdf(double x){
    return (1 + std::erf(x/std::pow(2, 0.5))) / 2;
}

double compute_d1(double S, double K, double T, double sigma, double rate){
    return (std::log(S/K) + ((rate + std::pow(sigma, 2) / 2) * T)) / (sigma * std::pow(T, 0.5));
}

double compute_d2(double d1, double T, double sigma){
    return d1 - sigma * std::pow(T, 0.5);
}

double compute_price(double S, double K, double T, double rate, double sigma){
    double d1 = compute_d1(S, K, T, sigma, rate);
    double d2 = compute_d2(d1, T, sigma);
    double N1 = normal_cdf(d1);
    double N2 = normal_cdf(d2);
    return S * N1 - K * std::exp((-rate * T)) * N2;
}

#endif