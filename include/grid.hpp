#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <cassert>

enum CellType { INTERIOR, FIXED, HOLE };

class Grid { 
private:
    int rows;
    int cols;
    std::vector<double> temps;
    std::vector<CellType> cellTypes;
public:
    Grid(int r, int c){
        assert(r > 0 && c > 0);
        rows = r;
        cols = c;
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temps.push_back(0);
                if (i == 0 || i == r - 1 || j == 0 || j == c - 1){cellTypes.push_back(FIXED);}
                else{cellTypes.push_back(INTERIOR);}
            }
        }
    }

    double& at(int r, int c){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return temps[r * cols + c];
    }

    const double& at(int r, int c) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return temps[r * cols + c];
    }

    int getRows() const {return rows;}
    int getCols() const {return cols;}

    void setType(int i, int j, CellType type){
        assert(i >= 0 && i < getRows() && j >= 0 && j < getCols());
        cellTypes[i * cols + j] = type;
    }

    CellType getType(int i, int j) const {
        assert(i >= 0 && i < getRows() && j >= 0 && j < getCols());
        return cellTypes[i * cols + j];
    }

    void maskRect(int i0, int j0, int i1, int j1, CellType type){
        for (int i = i0; i < i1; i++){
            for (int j = j0; j < j1; j++){
                setType(i, j, type);
            }
        }
    }

    double* getTempsPtr(){return temps.data();}
    CellType* getTypesPtr(){return cellTypes.data();}
};

#endif