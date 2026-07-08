#include <iostream>
#include <cmath>
#include "grid.hpp"
#include "sor.hpp"

const double PI = 3.14159265358979323846;

int main(){
    int sizes[] = {25, 50, 100};

    // Experiment 1: T(x,y) = x*y  (zero discretization error;
    // measures pure iteration error)
    std::cout << "T(x,y) = x*y\n";
    for (int s = 0; s < 3; s++){
        double h = 1.0/(sizes[s] - 1.0);
        Grid g(sizes[s], sizes[s]);
        for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){
                if (i == 0 || j == 0 || i == g.getRows() - 1 || j == g.getCols() - 1){
                    g.at(i,j) = (i * h) * (j * h);
                }
            }
        }
        sor_solve(g, 1e-8, 50000, 1.5);
        double maxError = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double error = std::abs(g.at(i,j) - (i * h) * (j * h));
                if (error > maxError){
                    maxError = error;
                }
            }
        }
        std::cout << "N: " << sizes[s] << " Max Error: " << maxError << "\n";
    }

    // Experiment 2: T(x,y) = sin(pi x) * sinh(pi y)  
    std::cout << "\nT(x,y) = sin(pi x) * sinh(pi y)\n";
    for (int s = 0; s < 3; s++){
        double h = 1.0/(sizes[s] - 1.0);
        Grid g(sizes[s], sizes[s]);
        for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){
                if (i == 0 || j == 0 || i == g.getRows() - 1 || j == g.getCols() - 1){
                    g.at(i,j) = std::sin(PI * j * h) * std::sinh(PI * i * h);
                }
            }
        }
        sor_solve(g, 1e-8, 50000, 1.5);
        double maxError = 0;
        for (int i = 1; i < g.getRows() - 1; i++){
            for (int j = 1; j < g.getCols() - 1; j++){
                double error = std::abs(g.at(i,j) - std::sin(PI * j * h) * std::sinh(PI * i * h));
                if (error > maxError){
                    maxError = error;
                }
            }
        }
        std::cout << "N: " << sizes[s] << " Max Error: " << maxError << "\n";
    }
    return 0;
}