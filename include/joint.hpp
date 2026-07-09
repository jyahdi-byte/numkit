#ifndef JOINT_HPP
#define JOINT_HPP

#include <vector>
#include <cmath>

class Joint{
private:
    double x,y;
    bool anchoredX, anchoredY;
public:
    Joint(double x0, double y0, bool anchX = false, bool anchY = false){
        x = x0;
        y = y0;
        anchoredX = anchX;
        anchoredY = anchY;
    }
    double getX() const {
        return x;
    }
    double getY() const {
        return y;
    }
    bool isAnchoredX() const {
        return anchoredX;
    }
    bool isAnchoredY() const {
        return anchoredY;
    }
};

#endif