#ifndef SPACE_TIME_LOG_HPP
#define SPACE_TIME_LOG_HPP

#include <vector>
#include <cassert>
#include "grid1d.hpp"
 
class SpaceTimeLog {
private:
    int rows;
    int cols;
    std::vector<double> values;
public:
    SpaceTimeLog(int r, int c){
        assert(r > 0 && c > 0);
        rows = r;
        cols = c;
        values.resize(rows * cols);
    }   

    void setRow(int r, const Grid1D& vals){
        assert(r >= 0 && r < rows && vals.getPoints() == cols);
        for (int i = 0; i < vals.getPoints(); i++){
            values[r * cols + i] = vals.at(i);
        }
    }

    const double& at(int r, int c) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return values[r * cols + c];
    }

    int getRows() const {return rows;}
    int getCols() const {return cols;}
};

#endif