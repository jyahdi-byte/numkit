#ifndef BAR_HPP
#define BAR_HPP

#include <vector>
#include <cmath>
#include "joint.hpp"

class Bar {
private:
    int jointA_index, jointB_index;
    double k;
public:
    Bar(int a_index, int b_index, double kConstant){
        jointA_index = a_index;
        jointB_index = b_index;
        k = kConstant;
    }

    double getXLength(const std::vector<Joint>& joints) const {
        return std::abs(joints[jointA_index].getX() - joints[jointB_index].getX());
    }

    double getYLength(const std::vector<Joint>& joints) const {
        return std::abs(joints[jointA_index].getY() - joints[jointB_index].getY());
    }

    double getLength(const std::vector<Joint>& joints) const {
       return std::sqrt(std::pow(getXLength(joints), 2) +std::pow(getYLength(joints),2));
    }

    double getTheta(const std::vector<Joint>& joints) const {
        double dx = joints[jointB_index].getX() - joints[jointA_index].getX();
        double dy = joints[jointB_index].getY() - joints[jointA_index].getY();
        return std::atan2(dy, dx);
    }
};

#endif