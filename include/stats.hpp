#ifndef STATS_HPP
#define STATS_HPP

#include <vector>
#include <cmath>

inline double mean(const std::vector<double>& vec){
    double sum  = 0;
    for (int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum / vec.size();
}

inline double var(const std::vector<double>& vec){
    double avg = mean(vec);
    double sum = 0;
    for (int i = 0; i < vec.size(); i++){
        sum += std::pow(vec[i] - avg, 2);
    }
    return sum / vec.size();
}

inline double dev(const std::vector<double>& vec){
    return std::sqrt(var(vec));
}

#endif