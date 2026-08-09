#ifndef GRID1D_HPP
#define GRID1D_HPP

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

class Grid1D {
private:
    int points;
    double dx;
    std::vector<double> values;
public:
    Grid1D(int p, double x){
        assert(p > 0);
        points = p;
        dx = std::abs(x);
        values.resize(p);
    }

    void operator = (const Grid1D& other){
        if (points != other.points){
            points = other.points;
            values.resize(points);  
        }  
        for (int i = 0; i < other.points; i++){
                at(i) = other.at(i);
        }
        dx = other.dx;
    }

    Grid1D(const Grid1D& other){
        points = other.points;
        dx = other.dx;
        values.resize(points);
        for (int i = 0; i < other.points; i++){
            at(i) = other.at(i);
        }
    }
 
    double& at(int x){
        // periodic wraparound: index -1 maps to points-1, index points maps to 0
        if (x < 0){return at(points + x);}
        else if (x >= points){return at(x - points);}
        return values[x];
    }

    const double& at(int x) const {
        // periodic wraparound: index -1 maps to points-1, index points maps to 0
        if (x < 0){return at(points + x);}
        else if (x >= points){return at(x - points);}
        return values[x];
    } 

    void print() const {
        for (int i = 0; i < points; i++){
            std::cout << at(i) << " ";
        }
        std::cout << "\n";
    }

    double sum() const {
        double sum = 0;
        for (int i = 0; i < points; i++){
            sum += at(i);
        }
        return sum;
    }

    int getPoints() const {return points;}
    double getDx() const {return dx;}
    double getLength() const {return points * dx;}
};

#endif