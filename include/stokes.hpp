#ifndef _STOKES_SOLVER_H_
#define _STOKES_SOLVER_H_

#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<fstream>
#include<sstream>
#include<set>
#include<map>
#include<chrono>
#include<algorithm>
#include "TextParser.h"
#include "pardiso_solver.h"

class STOKES_solver{
    double alpha, resistance, mu;
    std::vector<std::vector<int>> element_v, element_p;
    std::vector<std::vector<double>> node, pressure_node;
    std::map<int, int> pressure_element_transform;
    std::vector<double> phi;
    std::vector<double> u,v,pressure;
    std::vector<int> inlet_boundary_node;

    std::string base_input_dir;

    int numOfNodeInElmVelocity, numOfNodeInElmPressure;

    int numOfOmpThreads;

    std::vector<double> gauss_point, gauss_weight;

    PARDISO_solver PARDISO;
    public:
        void input_parameter_from_tp_file();
        void input_domain_info();
        void input_boundary_info();
        void redifine_pressure_node_element();
        void pressure_velocity_initialize();

    void gauss_point_setting();

    void main_stokes();
    void ShapeFunctionC2D8_dNdr(std::vector<std::vector<double>> &dNdr, double g1, double g2);
    void ShapeFunctionC2D8(std::vector<double> &N, double g1, double g2);
    void ShapeFunctionC2D4(std::vector<double> &N, double g1, double g2);
    void calc_dxdr(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr, std::vector<std::vector<int>> element, std::vector<std::vector<double>> x, int ic, int numOfNodeInElm);
    void calc_inverse_matrix2x2(std::vector<std::vector<double>> &drdx, std::vector<std::vector<double>> dxdr);
    void calc_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> drdx, std::vector<std::vector<double>> &dNdr, int numOfNodeInElm);
    double calc_determinant2x2(std::vector<std::vector<double>> A);

    void calc_diffusive_and_Darcy_matrix(int ic);
    void calc_pressure_matrix(int ic);

    void export_vtu_pressure(const std::string &file);
    void export_vtu_velocity(const std::string &file);




    public:
        TextParser tp;
};

#endif //_PARDISO_SOLVER_H_