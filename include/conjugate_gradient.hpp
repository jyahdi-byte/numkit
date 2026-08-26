#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include <cassert>
#include <vector>
#include <cmath>

#include "grid.hpp"

double dot(std::vector<double> a, std::vector<double> b){
    assert(a.size() == b.size());
    double result = 0;
    for (int i = 0; i < a.size(); i++){result += a[i] * b[i];}
    return result;
}

std::vector<double> compute_b(const Grid& g) {
    int rows = g.getRows();
    int cols = g.getCols();
    std::vector<double> b(rows * cols, 0.0);

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            if (g.getType(i, j) == INTERIOR) {
                double sum = 0;
                if (g.getType(i + 1, j) == FIXED) {sum += g.getFacesK(i, j, 0) * g.at(i + 1, j);}
                if (g.getType(i - 1, j) == FIXED) {sum += g.getFacesK(i, j, 1) * g.at(i - 1, j);}
                if (g.getType(i, j + 1) == FIXED) {sum += g.getFacesK(i, j, 2) * g.at(i, j + 1);}
                if (g.getType(i, j - 1) == FIXED) {sum += g.getFacesK(i, j, 3) * g.at(i, j - 1);}
                b[i * cols + j] = sum;
            }
        }
    }
    return b;
} 

std::vector<double> apply_A(const std::vector<double>& p, int rows, int cols, const Grid& A){
    std::vector<double> result(rows * cols, 0.0);
    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            if (A.getType(i, j) == INTERIOR) {
                std::vector<double> components = {p[(i+1) * cols + j], p[(i-1) * cols + j], p[i * cols + j + 1], p[i * cols + j - 1]}; 
                double sum = (A.getFacesK(i, j, 0) * components[0] + A.getFacesK(i, j, 1) * components[1] + A.getFacesK(i, j, 2) * components[2] + A.getFacesK(i, j, 3) * components[3]);
                double center = p[i * cols + j];
                result[i * cols + j] = A.getTotalK(i, j) * center - sum;
            }
        }
    }
    return result;
}

int cg_solve(Grid& A, double tol, int iterations){
    int sweeps = 0;
    std::vector<double> x(A.getRows() * A.getCols()); std::vector<double> r_new(A.getRows() * A.getCols());
    for (int i = 0; i < A.getRows(); i++)
        for (int j = 0; j < A.getCols(); j++){x[i * A.getCols() + j] = 0;}
    std::vector<double> b = compute_b(A);
    std::vector<double> r = b;
    std::vector<double> p = r;
    for (int k = 0; k < iterations; k++){
        std::vector<double> Ap = apply_A(p, A.getRows(), A.getCols(), A);
        double alpha = dot(r, r) / dot (p , Ap);
        for (int i = 0; i < x.size(); i++){x[i] += alpha * p[i];}
        sweeps++;
        for (int i = 0; i < r_new.size(); i++){r_new[i] = r[i] + (-alpha * Ap[i]);}
        double beta = dot(r_new, r_new) / dot(r, r);
        for (int i = 0; i < p.size(); i++){p[i] = r_new[i] + beta * p[i];}
        r = r_new; 
        if (tol > std::sqrt(dot(r,r))){
            for (int i = 0; i < A.getRows(); i++){
                for (int j = 0; j < A.getCols(); j++){
                    if (A.getType(i, j) == INTERIOR){
                        A.at(i, j) = x[i * A.getCols() + j];
                    }
                }
            }
            return sweeps;
        }
    }
    for (int i = 0; i < A.getRows(); i++){
        for (int j = 0; j < A.getCols(); j++){
            if (A.getType(i, j) == INTERIOR){
                A.at(i, j) = x[i * A.getCols() + j];
            }
        }
    }
    return sweeps;
}

int pcg_solve(Grid& A, double tol, int iterations){
    int sweeps = 0;
    std::vector<double> x(A.getRows() * A.getCols()); std::vector<double> r_new(A.getRows() * A.getCols()); 
    for (int i = 0; i < A.getRows(); i++)
        for (int j = 0; j < A.getCols(); j++){x[i * A.getCols() + j] = 0;}
    std::vector<double> b = compute_b(A);
    std::vector<double> r = b;
    std::vector<double> z = r;
    for (int i = 1; i < A.getRows() - 1; i++){
        for (int j = 1; j < A.getCols() - 1; j++){
            if (A.getType(i, j) == INTERIOR){
                z[i * A.getCols() + j] /= A.getTotalK(i, j);
            }
        } 
    }
    std::vector<double> p = z;
    for (int k = 0; k < iterations; k++){
        std::vector<double> Ap = apply_A(p, A.getRows(), A.getCols(), A);
        double alpha = dot(r, z) / dot (p , Ap);
        for (int i = 0; i < x.size(); i++){x[i] += alpha * p[i];}
        sweeps++;
        for (int i = 0; i < r_new.size(); i++){r_new[i] = r[i] + (-alpha * Ap[i]);}
        std::vector<double> z_new = r_new;
        for (int i = 1; i < A.getRows() - 1; i++){
            for (int j = 1; j < A.getCols() - 1; j++){
                if (A.getType(i, j) == INTERIOR){z_new[i * A.getCols() + j] /= A.getTotalK(i, j);}
            }
        }
        double beta = dot(r_new, z_new) / dot(r, z);
        for (int i = 0; i < p.size(); i++){p[i] = z_new[i] + beta * p[i];}
        r = r_new; z = z_new;
        if (tol > std::sqrt(dot(r,r))){
            for (int i = 0; i < A.getRows(); i++){
                for (int j = 0; j < A.getCols(); j++){
                    if (A.getType(i, j) == INTERIOR){A.at(i, j) = x[i * A.getCols() + j];}
                }
            }
            return sweeps;
        }
    }
    for (int i = 0; i < A.getRows(); i++){
        for (int j = 0; j < A.getCols(); j++){
            if (A.getType(i, j) == INTERIOR){
                A.at(i, j) = x[i * A.getCols() + j];
            }
        }
    }
    return sweeps;
}

#endif