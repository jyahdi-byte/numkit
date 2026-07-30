#include <iostream>
#include <cassert>

#include "grid.hpp"

int main(){
    int rows = 4;
    int cols = 4;
    Grid g(rows,cols);
    g.at(2,3) = 25.0;
    double x = g.at(2,3);
    assert(x == 25.0);
    assert(g.getRows() == rows);

    std::cout << "PASS\n";

    return 0;
}