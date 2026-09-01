#include <iostream>
#include <cassert>
#include <cmath>

#include "grid.hpp"
#include "slice_grid.hpp"

int main(){
    int rows = 6, cols = 7, depth = 8;
    Grid3D g(rows, cols, depth);

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            for (int k = 0; k < depth; k++){
                g.at(i,j,k) = (i + j) * k;
            }
        }
    }

    Grid sliceI = slice_grid3d(g, AXIS_I, 2);
    Grid sliceJ = slice_grid3d(g, AXIS_J, 3);
    Grid sliceK = slice_grid3d(g, AXIS_K, 4);

    assert(sliceI.getRows() == cols && sliceI.getCols() == depth);
    assert(sliceJ.getRows() == rows && sliceJ.getCols() == depth);
    assert(sliceK.getRows() == rows && sliceK.getCols() == cols);

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            assert(g.at(i, j, 4) == sliceK.at(i, j));
            assert(g.getType(i, j, 4) == sliceK.getType(i, j));
            assert(g.getConductivity(i, j, 4) == sliceK.getConductivity(i, j));
        }
    }

    for (int j = 0; j < cols; j++){
        for (int k = 0; k < depth; k++){
            assert(g.at(2, j, k) == sliceI.at(j, k));
            assert(g.getType(2, j, k) == sliceI.getType(j, k));
            assert(g.getConductivity(2, j, k) == sliceI.getConductivity(j, k));
        }
    }

    for (int i = 0; i < rows; i++){
        for (int k = 0; k < depth; k++){
            assert(g.at(i, 3, k) == sliceJ.at(i, k));
            assert(g.getType(i, 3, k) == sliceJ.getType(i, k));
            assert(g.getConductivity(i, 3, k) == sliceJ.getConductivity(i, k));
        }
    }

    std::cout << "PASS\n";

    return 0;
}