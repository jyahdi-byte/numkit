#include <iostream>
#include <vector>
#include <cmath>

#include "grid.hpp"
#include "gauss_seidel.hpp"
#include "sor.hpp"
#include "conjugate_gradient.hpp"

int main(){
    const char* solverNames[] = {"GS", "SOR", "CG"};
    std::vector<Grid> grids;
    for (int n = 10; n < 100; n *= 2){grids.push_back(Grid(n,n));}
    for (size_t i = 0; i < grids.size(); i++){
        for (int j = 0; j < grids[i].getCols(); j++){grids[i].at(0,j) = 10;}
    }

    std::vector<Grid> gs_grids;
    std::vector<Grid> sor_grids;
    std::vector<Grid> cg_grids;
    for (size_t i = 0; i < grids.size(); i++){
        gs_grids.push_back(grids[i]);
        sor_grids.push_back(grids[i]);
        cg_grids.push_back(grids[i]);
    }

    std::vector<int> gs_sweeps;
    std::vector<int> sor_sweeps;
    std::vector<int> cg_sweeps;
    for (size_t k = 0; k < grids.size(); k++){
        gs_sweeps.push_back(gauss_seidel_solve(gs_grids[k], 1e-10, 20000));
        sor_sweeps.push_back(sor_solve(sor_grids[k], 1e-10, 20000));
        cg_sweeps.push_back(cg_solve(cg_grids[k], 1e-10, 20000));
    }
    std::vector<std::vector<int>> all_sweeps = {gs_sweeps, sor_sweeps, cg_sweeps};
 
    std::cout << "SWEEP COUNTS\n";
    for (size_t i = 0; i < grids.size(); i++){
        std::cout << "\n" << "Grid Size: " << 10 * std::pow(2, i) << "\n";
        for (size_t j = 0; j < all_sweeps.size(); j++){
            std::cout << solverNames[j] << ": " << all_sweeps[j][i] << "\n";
        }
    }

    return 0;
}