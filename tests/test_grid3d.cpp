#include <iostream>
#include <cassert>
#include <cmath>

#include "grid.hpp"

int main(){
    int rows = 4;
    int cols = 4;
    int depth = 4;
    Grid3D g(rows,cols,depth);

    g.at(2,2,2) = 25.0;
    double x = g.at(2,2,2);
    assert(x == 25.0);
    assert(g.getRows() == rows);
    assert(g.getCols() == cols);
    assert(g.getDepth() == depth);

    assert(g.getType(0,2,2) == FIXED);
    assert(g.getType(0,0,0) == FIXED);
    assert(g.getType(1,1,1) == INTERIOR);

    g.maskRect(1,1,1,2,2,2, HOLE);
    assert(g.getType(1,1,1) == HOLE);

    assert(g.getType(0,3,3) == FIXED);
    assert(g.getType(2,2,2) == INTERIOR);

    for (int d = 0; d < depth; d++){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                if (g.getType(i,j,d) == INTERIOR){
                    double sum = 0;
                    for (int dir = 0; dir < 6; dir++){sum += g.getFacesK(i,j,d,dir);}
                    assert(std::abs(sum - g.getTotalK(i,j,d)) < 1e-12);
                }
            }
        }
    }

    std::cout << "PASS\n";

    return 0;
}