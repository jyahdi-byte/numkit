#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <cassert>

class Grid {
private:
    int rows;
    int cols;
    std::vector<double> temps;
public:
    Grid(int r, int c){
        assert(r > 0 && c > 0);
        rows = r;
        cols = c;
        for (int i = 0; i < r * c; i++){
            temps.push_back(0);
        }
    }

    double& at(int r, int c){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return temps[r * cols + c];
    }

    int getRows(){
        return rows;
    }
    int getCols(){
        return cols;
    }
};

#endif