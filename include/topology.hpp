#ifndef _TOPOLOGY_H_
#define _TOPOLOGY_H_

#include "stokes.hpp"
#include <nlopt.hpp>

class topology{
    public:
    STOKES_solver Stokes_main;
    STOKES_solver Stokes_adjoint;

    TextParser tp;

    double currentVolume, TotalVolume;

    int loop;

    std::vector<std::vector<double>> adjoint_force;
    std::vector<double> dCdr;
    std::vector<double> rho;
    std::vector<double> element_volume;

    void initialize();
    double return_objective_function(const std::vector<double> &rho_MMA, std::vector<double> &grad, void *my_func_data);
    
    double calc_C();
    void calc_sensivity();
    void calc_adjoint_force();


};

#endif