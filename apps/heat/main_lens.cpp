#include <iostream>
#include <cmath>
#include <thread>
#include <vector>

#include "grid.hpp"
#include "sor_rb_mt.hpp"
#include "ppm.hpp"

// sets conductivity k along a straight band from (r0,c0) to (r1,c1), given band half-thickness
void finBand(Grid& g, double r0, double c0, double r1, double c1, double halfThickness, double k){
    double dr = r1 - r0, dc = c1 - c0;
    double lineLenSq = dr * dr + dc * dc;
    for (int i = 1; i < g.getRows() - 1; i++){
        for (int j = 1; j < g.getCols() - 1; j++){
            if (g.getType(i, j) != INTERIOR) continue;
            double t = ((i - r0) * dr + (j - c0) * dc) / lineLenSq;
            if (t < 0 || t > 1) continue;
            double perpDist = std::abs((i - r0) * dc - (j - c0) * dr) / std::sqrt(lineLenSq);
            if (perpDist < halfThickness){
                g.setConductivity(i, j, k);
            }
        }
    }
}

// radial conductivity bump, center to edge
void bump(Grid& g, int cx, int cy, int r, double kCenter, double kEdge){
    for (int i = cx - r; i <= cx + r; i++){
        for (int j = cy - r; j <= cy + r; j++){
            if (i <= 0 || i >= g.getRows() - 1 || j <= 0 || j >= g.getCols() - 1) continue;
            if (g.getType(i, j) != INTERIOR) continue;
            double dist = std::sqrt(std::pow(i - cx, 2) + std::pow(j - cy, 2));
            if (dist <= r){
                double frac = dist / r;
                g.setConductivity(i, j, kCenter + frac * (kEdge - kCenter));
            }
        }
    }
}

int main(){
    int rows = 280, cols = 420;
    Grid g(rows, cols);

    for (int i = 0; i < rows; i++){
        g.setType(i, 0, HOLE);
        g.setType(i, cols - 1, HOLE);
    }
    for (int j = 0; j < cols; j++){
        g.setType(0, j, HOLE);
        g.setType(rows - 1, j, HOLE);
    }

    // low-conductivity background so heat is forced to travel mostly through the fins
    for (int i = 1; i < rows - 1; i++){
        for (int j = 1; j < cols - 1; j++){
            if (g.getType(i, j) == INTERIOR){
                g.setConductivity(i, j, 0.15);
            }
        }
    }

    // wide hot base, bottom-left
    for (int i = rows - 30; i < rows - 1; i++){
        for (int j = 1; j < 160; j++){
            g.setType(i, j, FIXED);
            g.at(i, j) = 100.0;
        }
    }

    // cold sink, top-right, catches the fin tips
    for (int i = 1; i < 25; i++){
        for (int j = 250; j < cols - 1; j++){
            g.setType(i, j, FIXED);
            g.at(i, j) = 0.0;
        }
    }

    // fan of full-conductivity fins from the hot base up to the cold sink
    int numFins = 7;
    for (int f = 0; f < numFins; f++){
        double c0 = 20 + f * 18;
        double c1 = 260 + f * 20;
        finBand(g, rows - 20, c0, 15, c1, 7.0, 1.0);
    }

    // a couple of conductivity bumps along the fins for extra visual variety
    bump(g, rows / 2, 140, 22, 1.0, 0.05);
    bump(g, rows / 2 - 40, 280, 18, 0.9, 1.0);

    write_ppm_conductivity(g, "lens_conductivity.ppm");
    write_ppm_hue(g, "lens_before.ppm");

    int num_threads = std::thread::hardware_concurrency();
    if (num_threads <= 0) num_threads = 4;
    int sweeps = sor_rb_mt_solve(g, 1e-6, 200000, num_threads);
    std::cout << sweeps << "\n";

    write_ppm_hue(g, "lens_after.ppm");
    write_ppm_composite(g, "lens_composite_after.ppm", 0.6); 

    return 0;
}