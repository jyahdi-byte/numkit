#ifndef PPM_HPP
#define PPM_HPP

#include "grid.hpp"
#include <string>
#include <fstream>
#include <cmath>

void write_ppm(const Grid& g, const std::string& filename){
    double tmax = g.at(0,0);
    double tmin = g.at(0,0);
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (g.at(i,j) > tmax){
                tmax = g.at(i,j);
            }
            else if (g.at(i,j) < tmin){
                tmin = g.at(i,j);
            }
        }
    }

    std::ofstream out(filename);
    out << "P3" << "\n" << g.getCols() << " " << g.getRows() << "\n" << "255" << "\n";
    double t;
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (tmax - tmin == 0){
                t = 0.5;
            }
            else{
                t = (g.at(i,j) - tmin)/(tmax - tmin);
            }
            out << std::lround(t * 255) << " 0 " << std::lround((1-t) * 255) << " ";
        }
        out << "\n";
    }
}


#endif