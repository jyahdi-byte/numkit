#include <iostream>
#include <vector>
#include <cmath>

#include "grid.hpp"
#include "gauss_seidel.hpp"
#include "sor.hpp"
#include "gauss_seidel_rb.hpp"
#include "sor_rb.hpp"

int main(){
    const char* solverNames[] = {"GS", "SOR", "GS_RB", "SOR_RB"};
    std::vector<Grid> grids;
    for (int n = 10; n < 100; n *= 2){grids.push_back(Grid(n,n));}
    for (size_t i = 0; i < grids.size(); i++){
        for (int j = 0; j < grids[i].getCols(); j++){grids[i].at(0,j) = 10;}
    }

    std::vector<Grid> gs_grids;
    std::vector<Grid> sor_grids;
    std::vector<Grid> gs_rb_grids;
    std::vector<Grid> sor_rb_grids;
    for (size_t i = 0; i < grids.size(); i++){
        gs_grids.push_back(grids[i]);
        sor_grids.push_back(grids[i]);
        gs_rb_grids.push_back(grids[i]);
        sor_rb_grids.push_back(grids[i]);
    }

    std::vector<int> gs_sweeps;
    std::vector<int> sor_sweeps;
    std::vector<int> gs_rb_sweeps;
    std::vector<int> sor_rb_sweeps;
    for (size_t k = 0; k < grids.size(); k++){
        gs_sweeps.push_back(gauss_seidel_solve(gs_grids[k], 1e-10, 20000));
        sor_sweeps.push_back(sor_solve(sor_grids[k], 1e-10, 20000));
        gs_rb_sweeps.push_back(gauss_seidel_rb_solve(gs_rb_grids[k], 1e-10, 20000));
        sor_rb_sweeps.push_back(sor_rb_solve(sor_rb_grids[k], 1e-10, 20000));
    }
    std::vector<std::vector<int>> all_sweeps = {gs_sweeps, sor_sweeps, gs_rb_sweeps, sor_rb_sweeps};

    std::cout << "SWEEP COUNTS\n";
    for (size_t i = 0; i < grids.size(); i++){
        std::cout << "\n" << "Grid Size: " << 10 * std::pow(2, i) << "\n";
        for (size_t j = 0; j < grids.size(); j++){
            std::cout << solverNames[j] << ": " << all_sweeps[j][i] << "\n";
        }
    }

    return 0;
}