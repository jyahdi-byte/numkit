#ifndef SLICE_GRID_HPP
#define SLICE_GRID_HPP

#include <cassert>

#include "grid.hpp"

enum SliceAxis { AXIS_I, AXIS_J, AXIS_K };

Grid slice_grid3d(const Grid3D& g, SliceAxis axis, int index){
    if (axis == AXIS_I){
        assert(index >= 0 && index < g.getRows());
        Grid g1(g.getCols(), g.getDepth());

        for (int j = 0; j < g.getCols(); j++){
            for (int k = 0; k < g.getDepth(); k++){
                g1.at(j,k) = g.at(index,j,k);
                g1.setType(j,k,g.getType(index,j,k));
                g1.setConductivity(j,k,g.getConductivity(index,j,k));
            }
        }
        return g1;
    }
    
    else if (axis == AXIS_J){
        assert(index >= 0 && index < g.getCols());
        Grid g1(g.getRows(), g.getDepth());

        for (int i = 0; i < g.getRows(); i++){
            for (int k = 0; k < g.getDepth(); k++){
                g1.at(i,k) = g.at(i,index,k);
                g1.setType(i,k,g.getType(i,index,k));
                g1.setConductivity(i,k,g.getConductivity(i,index,k));
            }
        }
        return g1;
    }

    else {
        assert(index >= 0 && index < g.getDepth());
        Grid g1(g.getRows(), g.getCols());

        for (int i = 0; i < g.getRows(); i++){
            for (int j = 0; j < g.getCols(); j++){
                g1.at(i,j) = g.at(i,j,index);
                g1.setType(i,j,g.getType(i,j,index));
                g1.setConductivity(i,j,g.getConductivity(i,j,index));
            }
        }
        return g1;
    }


}

#endif