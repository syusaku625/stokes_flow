#ifndef _MMA_H_
#define _MMA_H_

#include <iomanip> // std::setw(int), std::setfill(char)
#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include "topology.hpp"

class MMA{
    public:
        void simp_usingMMA(topology &simp);
        static double return_objective_function(const std::vector<double> &rho_MMA, std::vector<double> &grad, void *my_func_data);
        static double constraint_volume(const std::vector<double> &rho_MMA, std::vector<double> &grad, void *my_func_data);
};

#endif