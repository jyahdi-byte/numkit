#include <iostream>
#include <vector>
#include "joint.hpp"
#include "bar.hpp"

int main(){
    std::vector<Joint> joints;
    joints.push_back(Joint(0.0,0.0));
    joints.push_back(Joint(3.0,4.0));

    Bar bar(0, 1, 100);

    std::cout << "Length: " << bar.getLength(joints) << "\n";
    std::cout << "Theta: " << bar.getTheta(joints);

    return 0;
}