#include <iostream>
#include <cassert>
#include <cmath>
#include <numbers>
#include <string>

#include "grid.hpp"
#include "gauss_seidel_rb_mt.hpp"
#include "sor_rb_mt.hpp"
#include "gauss_seidel_rb.hpp"
#include "sor_rb.hpp"

enum SolverType { GAUSS_SEIDEL, SOR };
enum ExperimentType { EXP1, EXP2 };

double validate_rb_mt_solver(SolverType solver, int N, ExperimentType exp, double tolerance){
    double h = 1.0/(N - 1.0);
    Grid g(N, N);
    Grid g1(N, N);
    if (exp == EXP1){
        for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){
                if (i == 0 || j == 0 || i == g.getRows() - 1 || j == g.getCols() - 1){
                    g.at(i,j) = (i * h) * (j * h);
                }
            }
        }
    }
    else if (exp == EXP2){
        for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){
                if (i == 0 || j == 0 || i == g.getRows() - 1 || j == g.getCols() - 1){
                    g.at(i,j) = std::sin(std::numbers::pi * j * h) * std::sinh(std::numbers::pi * i * h);
                }
            }
        }
    }
    for (int i = N/4; i < 3*N/4; i++){
        for (int j = N/4; j < 3*N/4; j++){
            g.setType(i,j,HOLE);
        }
    }
    g1 = g;

    if (solver == GAUSS_SEIDEL){
        gauss_seidel_rb_mt_solve(g, 1e-10, 50000, 8);
        gauss_seidel_rb_solve(g1, 1e-10, 50000);
    }
    else if (solver == SOR){
        sor_rb_mt_solve(g, 1e-10, 50000, 8);
        sor_rb_solve(g1, 1e-10, 50000);
    }
    double maxError = 0;
    for (int i = 1; i < g.getRows() - 1; i++){
        for (int j = 1; j < g.getCols() - 1; j++){
            double error;
            if (g.getType(i,j) == HOLE || g.getType(i,j) == FIXED){continue;}
            error = std::abs(g.at(i,j) - g1.at(i,j));
            if (error > maxError){
                maxError = error;
            }
        }
    }
    return maxError;    
}

int main(){
    int sizes[] = {25, 50, 100};
    double tolerance = 1e-10;

    std::vector<std::string> solverNames = {"Gauss-Seidel Red-Black Multithreaded", "SOR Red-Black Multithreaded"};
    std::vector<std::string> experimentNames = {"T(x,y) = x*y", "T(x,y) = sin(pi x) * sinh(pi y)"};

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            std::cout << "\nValidating " << solverNames[i] << " for " << experimentNames[j] << "\n";
            for (int s = 0; s < 3; s++){
                double maxError = validate_rb_mt_solver((SolverType)i, sizes[s], (ExperimentType)j, tolerance);
                if (maxError > tolerance){
                    std::cout << "Error exceeds tolerance of " << tolerance << " for N = " << sizes[s] << "\n";
                    assert(false);
                }
            }
        }
    }

    return 0;
}