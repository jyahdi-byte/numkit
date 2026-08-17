#include <iostream>
#include <cassert>

#include "grid.hpp"

int main(){
    int rows = 4;
    int cols = 4;
    Grid g(rows,cols);

    // basic at() still works like before
    g.at(2,3) = 25.0;
    double x = g.at(2,3);
    assert(x == 25.0);
    assert(g.getRows() == rows);

    // outer ring should default to FIXED
    assert(g.getType(0, 2) == FIXED);

    // inside the ring should default to INTERIOR
    assert(g.getType(1, 2) == INTERIOR);

    // mask a single cell as a HOLE
    g.maskRect(1, 2, 2, 3, HOLE);
    assert(g.getType(1,2) == HOLE);

    // corner should still be FIXED, untouched by the mask
    assert(g.getType(0,3) == FIXED);

    // neighbor of the masked cell should NOT have been masked too
    assert(g.getType(2,2) == INTERIOR);

    std::cout << "PASS\n";

    return 0;
}