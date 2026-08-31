#include <iostream>
#include <cmath>
#include <thread>
#include <vector>

#include "grid.hpp"
#include "sor_rb_mt.hpp"
#include "slice_grid.hpp"
#include "ppm.hpp"

// sets conductivity k along a straight 3D band from (r0,c0,d0) to (r1,c1,d1), given band half-thickness
void finBand3d(Grid3D& g, double r0, double c0, double d0, double r1, double c1, double d1, double halfThickness, double k){
    double dr = r1 - r0, dc = c1 - c0, dd = d1 - d0;
    double lineLenSq = dr * dr + dc * dc + dd * dd;
    for (int i = 1; i < g.getRows() - 1; i++){
        for (int j = 1; j < g.getCols() - 1; j++){
            for (int m = 1; m < g.getDepth() - 1; m++){
                if (g.getType(i, j, m) != INTERIOR) continue;
                double t = ((i - r0) * dr + (j - c0) * dc + (m - d0) * dd) / lineLenSq;
                if (t < 0 || t > 1) continue;
                double pr = r0 + t * dr, pc = c0 + t * dc, pd = d0 + t * dd;
                double perpDist = std::sqrt(std::pow(i - pr, 2) + std::pow(j - pc, 2) + std::pow(m - pd, 2));
                if (perpDist < halfThickness){
                    g.setConductivity(i, j, m, k);
                }
            }
        }
    }
}

// radial (spherical) conductivity bump, center to edge
void bump3d(Grid3D& g, int cx, int cy, int cz, int r, double kCenter, double kEdge){
    for (int i = cx - r; i <= cx + r; i++){
        for (int j = cy - r; j <= cy + r; j++){
            for (int m = cz - r; m <= cz + r; m++){
                if (i <= 0 || i >= g.getRows() - 1 || j <= 0 || j >= g.getCols() - 1 || m <= 0 || m >= g.getDepth() - 1) continue;
                if (g.getType(i, j, m) != INTERIOR) continue;
                double dist = std::sqrt(std::pow(i - cx, 2) + std::pow(j - cy, 2) + std::pow(m - cz, 2));
                if (dist <= r){
                    double frac = dist / r;
                    g.setConductivity(i, j, m, kCenter + frac * (kEdge - kCenter));
                }
            }
        }
    }
}

int main(){
    int rows = 60, cols = 90, depth = 60;
    Grid3D g(rows, cols, depth);

    // outer shell is HOLE (insulated boundary), matching main_lens.cpp's border treatment
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            g.setType(i, j, 0, HOLE);
            g.setType(i, j, depth - 1, HOLE);
        }
    }
    for (int i = 0; i < rows; i++){
        for (int m = 0; m < depth; m++){
            g.setType(i, 0, m, HOLE);
            g.setType(i, cols - 1, m, HOLE);
        }
    }
    for (int j = 0; j < cols; j++){
        for (int m = 0; m < depth; m++){
            g.setType(0, j, m, HOLE);
            g.setType(rows - 1, j, m, HOLE);
        }
    }

    // low-conductivity background so heat is forced to travel mostly through the fins
    for (int i = 1; i < rows - 1; i++){
        for (int j = 1; j < cols - 1; j++){
            for (int m = 1; m < depth - 1; m++){
                if (g.getType(i, j, m) == INTERIOR){
                    g.setConductivity(i, j, m, 0.15);
                }
            }
        }
    }

    // wide hot base, bottom-left slab
    for (int i = rows - 12; i < rows - 1; i++){
        for (int j = 1; j < 35; j++){
            for (int m = 1; m < depth - 1; m++){
                g.setType(i, j, m, FIXED);
                g.at(i, j, m) = 100.0;
            }
        }
    }

    // cold sink, top-right slab, catches the fin tips
    for (int i = 1; i < 10; i++){
        for (int j = 55; j < cols - 1; j++){
            for (int m = 1; m < depth - 1; m++){
                g.setType(i, j, m, FIXED);
                g.at(i, j, m) = 0.0;
            }
        }
    }

    // fan of full-conductivity fins from the hot base up to the cold sink,
    // spread across depth so the fan isn't flat
    int numFins = 5;
    for (int f = 0; f < numFins; f++){
        double c0 = 5 + f * 8;
        double c1 = 60 + f * 6;
        double d0 = 10 + f * 8;
        double d1 = depth - 10 - f * 6;
        finBand3d(g, rows - 5, c0, d0, 4, c1, d1, 3.0, 1.0);
    }

    // conductivity bumps for extra visual variety
    bump3d(g, rows / 2, 30, depth / 2, 8, 1.0, 0.05);
    bump3d(g, rows / 2 - 8, 60, depth / 2 + 5, 6, 0.9, 1.0);

    // slice before solving: constant-k slice through the mid-depth, "looking down z"
    write_ppm_conductivity(slice_grid3d(g, AXIS_K, depth / 2), "lens3d_conductivity_k.ppm");
    write_ppm_hue(slice_grid3d(g, AXIS_K, depth / 2), "lens3d_before_k.ppm");

    int num_threads = std::thread::hardware_concurrency();
    if (num_threads <= 0) num_threads = 4;
    int sweeps = sor_rb_mt_solve(g, 1e-6, 200000, num_threads);
    std::cout << sweeps << "\n";

    // after solving: dump slices along all three axes to sanity-check the field
    write_ppm_hue(slice_grid3d(g, AXIS_K, depth / 4), "lens3d_after_k_q1.ppm");
    write_ppm_hue(slice_grid3d(g, AXIS_K, depth / 2), "lens3d_after_k_mid.ppm");
    write_ppm_hue(slice_grid3d(g, AXIS_K, 3 * depth / 4), "lens3d_after_k_q3.ppm");

    write_ppm_hue(slice_grid3d(g, AXIS_I, rows / 2), "lens3d_after_i_mid.ppm");
    write_ppm_hue(slice_grid3d(g, AXIS_J, cols / 2), "lens3d_after_j_mid.ppm");

    write_ppm_composite(slice_grid3d(g, AXIS_K, depth / 2), "lens3d_composite_after_k.ppm", 0.6);

    return 0;
}