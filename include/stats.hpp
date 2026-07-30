#ifndef STATS_HPP
#define STATS_HPP

#include <vector>
#include <cmath>
#include <cassert>

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

double RMSE(const std::vector<double>& points, const std::vector<double>& exacts){
    assert(points.size() == exacts.size());
    double sum  = 0;
    for (int i = 0; i < points.size(); i++){sum += std::pow(points[i] - exacts[i], 2);}
    return std::pow(sum/points.size(), 0.5);
}

double max_index(std::vector<double>& vec){
    double max = vec[0];
    double result = 0;
    for (int i = 0; i < vec.size(); i++){
        if (vec[i] > max){
            max = vec[i];
            result = i;
        }
    }
    return result;
}

#endif