#ifndef GRID_HPP
#define GRID_HPP

#include <vector> 
#include <cassert>
#include <cmath>

enum CellType { INTERIOR, FIXED, HOLE };

class Grid { 
private:
    int rows;
    int cols;
    std::vector<double> temps;
    std::vector<CellType> cellTypes;

    std::vector<double> conductivities;
    std::vector<double> faces_k; // flat: idx*4 + {0=up,1=down,2=right,3=left}
    std::vector<double> total_ks;
    double epsilon = 1e-9;
    std::vector<unsigned char> active; // 0 = false, 1 = true (vector<bool> has no .data())

    // --- conductivity-field internals ---

    double k_sum(int r1, int c1, int r2, int c2){
        assert(r1 < rows && r1 >= 0 && c1 < cols && c1 >= 0);
        assert(r2 < rows && r2 >= 0 && c2 < cols && c2 >= 0);
        return conductivities[r1 * cols + c1] + conductivities[r2 * cols + c2];
    }

    void update_k_component(int r, int c, int dr, int dc){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        assert(dr == 0 || dc == 0);
        
        if (dr == 1 && r + 1 < rows){faces_k[((r + 1) * cols + c) * 4 + 1] = faces_k[(r * cols + c) * 4 + 0];}
        else if (dr == -1 && r - 1 >= 0){faces_k[((r - 1) * cols + c) * 4 + 0] = faces_k[(r * cols + c) * 4 + 1];}
        else if (dc == 1 && c + 1 < cols){faces_k[(r * cols + (c + 1)) * 4 + 3] = faces_k[(r * cols + c) * 4 + 2];}
        else if (dc == -1 && c - 1 >= 0){faces_k[(r * cols + (c - 1)) * 4 + 2] = faces_k[(r * cols + c) * 4 + 3];}
    }

    void initialize_k(){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                double k_up; double k_down; double k_left; double k_right;
                if (i + 1 == rows ){k_up = conductivities[i * cols + j];}
                else if (k_sum(i, j, i + 1, j) == 0){k_up = 0;} else {k_up = (2 * conductivities[i * cols + j] * conductivities[(i + 1) * cols + j]) / k_sum(i, j, i + 1, j);}
                if (i - 1 < 0){k_down = conductivities[i * cols + j];}
                else if (k_sum(i, j, i - 1, j) == 0){k_down = 0;} else {k_down = (2 * conductivities[i * cols + j] * conductivities[(i - 1) * cols + j]) / k_sum(i, j, i - 1, j);}
                if (j + 1 == cols){k_right = conductivities[i * cols + j];}
                else if (k_sum(i, j, i, j + 1) == 0){k_right = 0;} else {k_right = (2 * conductivities[i * cols + j] * conductivities[i * cols + (j + 1)]) / k_sum(i, j, i, j + 1);}
                if (j - 1 < 0){k_left = conductivities[i * cols + j];}
                else if (k_sum(i, j, i, j - 1) == 0){k_left = 0;} else {k_left = (2 * conductivities[i * cols + j] * conductivities[i * cols + (j - 1)]) / k_sum(i, j, i, j - 1);}

                faces_k.push_back(k_up);
                faces_k.push_back(k_down);
                faces_k.push_back(k_right);
                faces_k.push_back(k_left);
                total_ks.push_back(k_up + k_down + k_right + k_left);
            }
        }
    }

    void initialize_active(){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                if (total_ks[i * cols + j] > epsilon){active.push_back(1);}
                else{active.push_back(0);}
            }
        }
    }

    void update_active(int r, int c){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        if (total_ks[r * cols + c] > epsilon){active[r * cols + c] = 1;}
        else{active[r * cols + c] = 0;}
    }

    void update_k(int r, int c){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        double k_up; double k_down; double k_left; double k_right;
        if (r + 1 == rows ){k_up = conductivities[r * cols + c];}
        else if (k_sum(r, c, r + 1, c) == 0){k_up = 0;} else {k_up = (2 * conductivities[r * cols + c] * conductivities[(r + 1) * cols + c]) / k_sum(r, c, r + 1, c);}
        if (r - 1 < 0){k_down = conductivities[r * cols + c];}
        else if (k_sum(r, c, r - 1, c) == 0){k_down = 0;} else {k_down = (2 * conductivities[r * cols + c] * conductivities[(r - 1) * cols + c]) / k_sum(r, c, r - 1, c);}
        if (c + 1 == cols){k_right = conductivities[r * cols + c];}
        else if (k_sum(r, c, r, c + 1) == 0){k_right = 0;} else {k_right = (2 * conductivities[r * cols + c] * conductivities[r * cols + (c + 1)]) / k_sum(r, c, r, c + 1);}
        if (c - 1 < 0){k_left = conductivities[r * cols + c];}
        else if (k_sum(r, c, r, c - 1) == 0){k_left = 0;} else {k_left = (2 * conductivities[r * cols + c] * conductivities[r * cols + (c - 1)]) / k_sum(r, c, r, c - 1);}

        faces_k[(r * cols + c) * 4 + 0] = k_up;
        faces_k[(r * cols + c) * 4 + 1] = k_down;
        faces_k[(r * cols + c) * 4 + 2] = k_right;
        faces_k[(r * cols + c) * 4 + 3] = k_left;

        update_k_component(r, c, 1, 0);
        update_k_component(r, c, -1, 0);
        update_k_component(r, c, 0, 1);
        update_k_component(r, c, 0, -1);

        total_ks[r * cols + c] = k_up + k_down + k_right + k_left;
        update_active(r, c);

        if (r + 1 < rows){
            total_ks[(r + 1) * cols + c] = faces_k[((r + 1) * cols + c) * 4 + 0] + faces_k[((r + 1) * cols + c) * 4 + 1] + faces_k[((r + 1) * cols + c) * 4 + 2] + faces_k[((r + 1) * cols + c) * 4 + 3];
            update_active(r + 1, c);
        }
        if (r - 1 >= 0){
            total_ks[(r - 1) * cols + c] = faces_k[((r - 1) * cols + c) * 4 + 0] + faces_k[((r - 1) * cols + c) * 4 + 1] + faces_k[((r - 1) * cols + c) * 4 + 2] + faces_k[((r - 1) * cols + c) * 4 + 3];
            update_active(r - 1, c);
        }
        if (c + 1 < cols){
            total_ks[r * cols + (c + 1)] = faces_k[(r * cols + (c + 1)) * 4 + 0] + faces_k[(r * cols + (c + 1)) * 4 + 1] + faces_k[(r * cols + (c + 1)) * 4 + 2] + faces_k[(r * cols + (c + 1)) * 4 + 3];
            update_active(r, c + 1);
        }
        if (c - 1 >= 0){
            total_ks[r * cols + (c - 1)] = faces_k[(r * cols + (c - 1)) * 4 + 0] + faces_k[(r * cols + (c - 1)) * 4 + 1] + faces_k[(r * cols + (c - 1)) * 4 + 2] + faces_k[(r * cols + (c - 1)) * 4 + 3];
            update_active(r, c - 1);
        }

    }

public:

    // --- construction ---

    Grid(int r, int c){
        assert(r > 0 && c > 0);
        rows = r;
        cols = c;
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                temps.push_back(0);
                if (i == 0 || i == r - 1 || j == 0 || j == c - 1){cellTypes.push_back(FIXED);}
                else{cellTypes.push_back(INTERIOR);}
                conductivities.push_back(1); 
            }
        }
        initialize_k(); initialize_active();
    }

    // --- grid shape ---

    int getRows() const {return rows;}
    int getCols() const {return cols;}

    // --- temperature access ---

    double& at(int r, int c){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return temps[r * cols + c];
    }

    const double& at(int r, int c) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return temps[r * cols + c];
    }

    double* getTempsPtr(){return temps.data();}
    const double* getTempsPtr() const {return temps.data();}

    // --- cell type ---

    void setType(int i, int j, CellType type){
        assert(i >= 0 && i < getRows() && j >= 0 && j < getCols());
        cellTypes[i * cols + j] = type;
        if (type == HOLE){
            conductivities[i * cols + j] = 0;
            update_k(i, j);
        }
        else if (type == FIXED){
            conductivities[i * cols + j] = 1;
            update_k(i, j);
        }
        update_active(i, j);
    } 

    CellType getType(int i, int j) const {
        assert(i >= 0 && i < getRows() && j >= 0 && j < getCols());
        return cellTypes[i * cols + j];
    }

    CellType* getTypesPtr(){return cellTypes.data();}

    // --- conductivity field access ---

    void setConductivity(int r, int c, double value){
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        conductivities[r * cols + c] = value;
        update_k(r, c);
    }

    const double& getConductivity(int r, int c) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return conductivities[r * cols + c];
    }

    double getFacesK(int r, int c, int dir) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        assert(dir >= 0 && dir < 4);
        return faces_k[(r * cols + c) * 4 + dir];
    }

    double getTotalK(int r, int c) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return total_ks[r * cols + c];
    }

    bool getActive(int r, int c) const {
        assert(r < rows && r >= 0 && c < cols && c >= 0);
        return active[r * cols + c] != 0;
    }

    const double* getFacesKPtr() const {return faces_k.data();}
    const double* getTotalKPtr() const {return total_ks.data();}
    const unsigned char* getActivePtr() const {return active.data();}

    // --- masking ---

    void maskRect(int i0, int j0, int i1, int j1, CellType type){
        for (int i = i0; i < i1; i++){
            for (int j = j0; j < j1; j++){
                if (i >= 0 && i < rows && j >= 0 && j < cols){setType(i, j, type);}
            }
        }
    }

    void maskCircle(int ci, int cj, int r, CellType type){
        for (int i = ci - r; i < ci + r; i++){
            for (int j = cj - r; j < cj + r; j++){
                if (i >= 0 && i < rows && j >= 0 && j < cols){
                    if (std::pow(i - ci, 2) + std::pow(j - cj, 2) <= std::pow(r, 2)){setType(i, j, type);}
                }
            }
        }
    }
};

#endif