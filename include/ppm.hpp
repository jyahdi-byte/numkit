#ifndef PPM_HPP
#define PPM_HPP

#include <string>
#include <fstream> 
#include <cmath>
#include <vector>

#include "grid.hpp"

std::vector<int> hsv_to_rgb(double h, double s, double v){
    double c = v * s;
    double x = c * (1 - std::abs((std::fmod(h / 60.0, 2)) - 1));
    double m = v - c;
    std::vector<double> rgb(3);

    if (h >= 0 && h < 60){rgb = {c, x, 0};}
    else if (h >= 60 && h < 120){rgb = {x, c, 0};}
    else if (h >= 120 && h < 180){rgb = {0, c, x};}
    else if (h >= 180 && h < 240){rgb = {0, x, c};}
    else if (h >= 240 && h < 300){rgb = {x, 0, c};}
    else {rgb = {c, 0, x};}

    return {std::lround(255 * (m + rgb[0])), std::lround(255 * (m + rgb[1])), std::lround(255 * (m + rgb[2]))};
}

void write_ppm_hue(const Grid& g, const std::string& filename){
    double tmax, tmin;
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (g.getActive(i,j)){
                tmax = g.at(i,j); tmin = g.at(i,j);
                break;
            }
        }
    }
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (g.getActive(i, j) != 0){
                if (g.at(i,j) > tmax){tmax = g.at(i,j);}
                else if (g.at(i,j) < tmin){tmin = g.at(i,j);}
            }
        }
    }

    std::ofstream out(filename);
    out << "P3" << "\n" << g.getCols() << " " << g.getRows() << "\n" << "255" << "\n";
    double h, t;
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (g.getActive(i,j) == 0){out << "0 0 0 ";}
            else{
                if (tmax - tmin == 0){t = 0.5;}
                else{t = (g.at(i,j) - tmin)/(tmax - tmin);}
                h = 240 * (1 - t);
                double s = 0.95;
                double v = 0.95;
                std::vector<int> rgb = hsv_to_rgb(h, s, v);
                out << rgb[0] << " " << rgb[1] << " " << rgb[2] << " ";
            }
        }
        out << "\n";
    }
}

void write_ppm_conductivity(const Grid& g, const std::string& filename){
    std::ofstream out(filename);
    out << "P3" << "\n" << g.getCols() << " " << g.getRows() << "\n" << "255" << "\n";
    double c;
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            c = g.getConductivity(i,j);
            out << std::lround(c * 255) << " " << std::lround(c * 255) << " " << std::lround(c * 255) << " ";
        }
        out << "\n";
    }
}

void write_ppm_composite(const Grid& g, const std::string& filename, double gamma){
    double tmax, tmin;
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (g.getActive(i,j)){
                tmax = g.at(i,j); tmin = g.at(i,j);
                break;
            }
        }
    }
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (g.getActive(i, j) != 0){
                if (g.at(i,j) > tmax){tmax = g.at(i,j);}
                else if (g.at(i,j) < tmin){tmin = g.at(i,j);}
            }
        }
    }

    std::ofstream out(filename);
    out << "P3" << "\n" << g.getCols() << " " << g.getRows() << "\n" << "255" << "\n";
    double h, s, v, t;
    for (int i = 0; i < g.getRows(); i++){
        for (int j = 0; j < g.getCols(); j++){
            if (tmax - tmin == 0){t = 0.5;}
            else{t = (g.at(i,j) - tmin)/(tmax - tmin);}
            h = 240 * (1 - t);
            s = 0.95 * std::pow(g.getConductivity(i,j), gamma);
            v = 0.95 * std::pow(g.getConductivity(i, j), gamma);
            std::vector<int> rgb = hsv_to_rgb(h, s, v);
            out << rgb[0] << " " << rgb[1] << " " << rgb[2] << " ";
        }
        out << "\n";
    }
}

#endif