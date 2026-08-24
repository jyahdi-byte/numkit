#ifndef PPM_HPP
#define PPM_HPP

#include "grid.hpp"
#include <string>
#include <fstream> 
#include <cmath>

template <typename T>

void write_ppm(const T& data, const std::string& filename){
    double tmax = data.at(0,0);
    double tmin = data.at(0,0);
    for (int i = 0; i < data.getRows(); i++){
        for (int j = 0; j < data.getCols(); j++){
            if (data.at(i,j) > tmax){
                tmax = data.at(i,j);
            }
            else if (data.at(i,j) < tmin){
                tmin = data.at(i,j);
            }
        }
    }

    std::ofstream out(filename);
    out << "P3" << "\n" << data.getCols() << " " << data.getRows() << "\n" << "255" << "\n";
    double t;
    for (int i = 0; i < data.getRows(); i++){
        for (int j = 0; j < data.getCols(); j++){
            if (tmax - tmin == 0){
                t = 0.5;
            }
            else{
                t = (data.at(i,j) - tmin)/(tmax - tmin);
            }
            out << std::lround(t * 255) << " 0 " << std::lround((1-t) * 255) << " ";
        }
        out << "\n";
    }
}

#endif