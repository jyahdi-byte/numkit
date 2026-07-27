#ifndef GRID1D_HPP
#define GRID1D_HPP

#include <vector>
#include <cassert>
#include <cmath>

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

    int getPoints() const {return points;}
    double getDx() const {return dx;}
    double getLength() const {return points * dx;}
};

#endif