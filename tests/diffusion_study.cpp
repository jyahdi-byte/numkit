#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "grid1d.hpp"
#include "advection.hpp"
#include "space_time_log.hpp"
#include "stats.hpp"
#include "ppm.hpp"

int main(){
    int iterations = 20;
    Grid1D u(50, 0.02);
    std::vector<double> uvec;
    for (int i = 0; i < u.getPoints(); i++){u.at(i) = 0;}
    for (int i = 5; i < 9; i ++){u.at(i) = 1;}

    Grid1D v = u;
    SpaceTimeLog log(iterations + 1, u.getPoints());
    advection_solve(v, iterations, 1, 1, log);
    for (int i = 0; i < u.getPoints(); i++){uvec.push_back(v.at(i));}
    

    std::vector<std::string> filenames;
    std::vector<SpaceTimeLog> logs;
    for (int k = 0; k < 5; k ++){
        double cn = 0.1 + 0.2 * k;
        filenames.push_back("advecton_map");
        filenames[k] += std::to_string(cn);
        filenames[k] +=  ".ppm";
        logs.push_back(SpaceTimeLog(iterations + 1, u.getPoints()));
    }
    
    std::vector<double> RMSEs;
    for (int k = 0; k < filenames.size(); k++){
        Grid1D u0 = u;
        double cn = 0.1 + 0.2 * k;
        advection_solve(u0, iterations, 1, cn, logs[k]);
        std::vector<double> results;
        for (int j = 0; j < u.getPoints(); j++){results.push_back(u0.at(j));}
        RMSEs.push_back(RMSE(results, uvec));
        write_ppm(logs[k], filenames[k]);
    }


    std::vector<double> cns = {0.1, 0.3, 0.5, 0.7, 0.9};
    for (int k = 0; k < filenames.size(); k++){
        int idx = max_index(RMSEs);
        std::cout << "Courant Number: " << cns[idx] << " RMSE " << RMSEs[idx] << "\n";
        cns.erase(cns.begin() + idx);
        RMSEs.erase(RMSEs.begin() + idx);
    }

    return 0;
}