#include <iostream>
#include "grid.hpp"

int main(){
    Grid g(4,4);
    g.at(2,3) = 25.0;
    double x = g.at(2,3);
    std::cout << x << "\n";
    std::cout << g.getRows() << "\n";
    g.at(500,3);

    return 0;
}