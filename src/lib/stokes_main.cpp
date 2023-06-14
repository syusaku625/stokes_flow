#include"stokes.hpp"

using namespace std;

void STOKES_solver::ShapeFunctionC2D8_dNdr(vector<vector<double>> &dNdr, double g1, double g2)
{
    dNdr[0][0] = 2.5e-1 * (1e+0-g2) * (2e+0*g1 +      g2);
    dNdr[0][1] = 2.5e-1 * (1e+0-g1) * (     g1 + 2e+0*g2);
    dNdr[1][0] = 2.5e-1 * (1e+0-g2) * (2e+0*g1 -      g2);
    dNdr[1][1] = 2.5e-1 * (1e+0+g1) * (    -g1 + 2e+0*g2);
    dNdr[2][0] = 2.5e-1 * (1e+0+g2) * (2e+0*g1 +      g2);
    dNdr[2][1] = 2.5e-1 * (1e+0+g1) * (     g1 + 2e+0*g2);
    dNdr[3][0] = 2.5e-1 * (1e+0+g2) * (2e+0*g1 -      g2);
    dNdr[3][1] = 2.5e-1 * (1e+0-g1) * (    -g1 + 2e+0*g2);
    dNdr[4][0] = -g1 * (1e+0-g2);
    dNdr[4][1] = -5e-1 * (1e+0-g1*g1);
    dNdr[5][0] = 5e-1 * (1e+0-g2*g2);
    dNdr[5][1] = -g2 * (1e+0+g1);
    dNdr[6][0] = -g1 * (1e+0+g2);
    dNdr[6][1] = 5e-1 * (1e+0-g1*g1);
    dNdr[7][0] = -5e-1 * (1e+0-g2*g2);
    dNdr[7][1] = -g2 * (1e+0-g1);
}

void STOKES_solver::ShapeFunctionC2D8(std::vector<double> &N, double g1, double g2)
{
    N[0] = 2.5e-1 * (1e+0-g1) * (1e+0-g2) * (-1e+0-g1-g2);
    N[1] = 2.5e-1 * (1e+0+g1) * (1e+0-g2) * (-1e+0+g1-g2);
    N[2] = 2.5e-1 * (1e+0+g1) * (1e+0+g2) * (-1e+0+g1+g2);
    N[3] = 2.5e-1 * (1e+0-g1) * (1e+0+g2) * (-1e+0-g1+g2);
    N[4] = 5e-1   * (1e+0-g1*g1) * (1e+0-g2);
    N[5] = 5e-1   * (1e+0+g1)    * (1e+0-g2*g2);
    N[6] = 5e-1   * (1e+0-g1*g1) * (1e+0+g2);
    N[7] = 5e-1   * (1e+0-g1)    * (1e+0-g2*g2);
}

void STOKES_solver::ShapeFunctionC2D4(std::vector<double> &N, double g1, double g2)
{
    N[0] = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
    N[1] = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
    N[2] = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
    N[3] = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
}

void STOKES_solver::calc_dxdr(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr, vector<vector<int>> element, vector<vector<double>> x, int ic, int numOfNodeInElm)
{
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            dxdr[i][j] = 0e0;
            for (int k = 0; k < numOfNodeInElm; k++){
                dxdr[i][j] += dNdr[k][j] * x[element[ic][k]][i];
            }
        }
    }
}

void STOKES_solver::calc_inverse_matrix2x2(std::vector<std::vector<double>> &drdx, std::vector<std::vector<double>> dxdr)
{
    double inverse=calc_determinant2x2(dxdr);
    drdx[0][0]=1e0/inverse*dxdr[1][1];
    drdx[1][1]=1e0/inverse*dxdr[0][0];
    drdx[0][1]=1e0/inverse*(dxdr[0][1])*(-1);
    drdx[1][0]=1e0/inverse*(dxdr[1][0])*(-1);
}

void STOKES_solver::calc_dNdx(vector<vector<double>> &dNdx, vector<vector<double>> drdx, vector<vector<double>> &dNdr, int numOfNodeInElm)
{
    for (int i = 0; i < numOfNodeInElm; i++){
        for (int j = 0; j < 2; j++){
            dNdx[i][j] = 0e0;
            for (int k = 0; k < 2; k++){
                dNdx[i][j] += dNdr[i][k] * drdx[k][j];
            }
        }
    }
}

double STOKES_solver::calc_determinant2x2(std::vector<std::vector<double>> A)
{
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];
}

void STOKES_solver::calc_diffusive_and_Darcy_matrix(int ic)
{
    vector<vector<double>>  dNdr(numOfNodeInElmVelocity, vector<double>(2, 0.0));
    vector<double> Np(numOfNodeInElmPressure, 0.0), Nv(numOfNodeInElmVelocity, 0.0);
    vector<vector<double>> dxdr(2, vector<double>(2, 0.0));
    vector<vector<double>> drdx(2, vector<double>(2, 0.0));
    vector<vector<double>> dNdx(numOfNodeInElmVelocity, vector<double>(2, 0.0));
    for(int j=0; j<gauss_point.size(); j++){
        for(int k=0; k<gauss_point.size(); k++){
            ShapeFunctionC2D8_dNdr(dNdr, gauss_point[j], gauss_point[k]);
            ShapeFunctionC2D8(Nv, gauss_point[j], gauss_point[k]);
            ShapeFunctionC2D4(Np, gauss_point[j], gauss_point[k]);
            calc_dxdr(dxdr, dNdr, element_v, node, ic, numOfNodeInElmVelocity);
            calc_inverse_matrix2x2(drdx,dxdr);
            calc_dNdx(dNdx, drdx, dNdr, numOfNodeInElmVelocity);
            double partial_volume = calc_determinant2x2(dxdr);
            //6x6 matrix
            for(int l=0; l<numOfNodeInElmVelocity; l++){
                for(int m=0; m<numOfNodeInElmVelocity; m++){
                    for(int n=0; n<2; n++){
                        PARDISO.coo_add(make_pair(element_v[ic][l], element_v[ic][m]), -mu * dNdx[l][n]*dNdx[m][n] * partial_volume * gauss_weight[j] * gauss_weight[k]);
                        PARDISO.coo_add(make_pair(element_v[ic][l]+node.size(), element_v[ic][m]+node.size()), -mu * dNdx[l][n]*dNdx[m][n] * partial_volume * gauss_weight[j] * gauss_weight[k]);
                        PARDISO.coo_add(make_pair(element_v[ic][l], element_v[ic][m]), -resistance * alpha * (1e0 - phi[ic]) / (alpha + phi[ic]) * Nv[l] * Nv[m] * partial_volume * gauss_weight[j] * gauss_weight[k]);
                        PARDISO.coo_add(make_pair(element_v[ic][l]+node.size(), element_v[ic][m]+node.size()), -resistance * alpha * (1e0 - phi[ic]) / (alpha + phi[ic]) * Nv[l] * Nv[m] * partial_volume * gauss_weight[j] * gauss_weight[k]);
                    }
                }
            }
        }
    }
}

void STOKES_solver::calc_pressure_matrix(int ic)
{
    vector<vector<double>>  dNdr(numOfNodeInElmVelocity, vector<double>(2, 0.0));
    vector<double> Np(numOfNodeInElmPressure, 0.0), Nv(numOfNodeInElmVelocity, 0.0);
    vector<vector<double>> dxdr(2, vector<double>(2, 0.0));
    vector<vector<double>> drdx(2, vector<double>(2, 0.0));
    vector<vector<double>> dNdx(numOfNodeInElmVelocity, vector<double>(2, 0.0));
    for(int j=0; j<gauss_point.size(); j++){
        for(int k=0; k<gauss_point.size(); k++){
            ShapeFunctionC2D8_dNdr(dNdr, gauss_point[j], gauss_point[k]);
            ShapeFunctionC2D8(Nv, gauss_point[j], gauss_point[k]);
            ShapeFunctionC2D4(Np, gauss_point[j], gauss_point[k]);
            calc_dxdr(dxdr, dNdr, element_v, node, ic, numOfNodeInElmVelocity);
            calc_inverse_matrix2x2(drdx,dxdr);
            calc_dNdx(dNdx, drdx, dNdr, numOfNodeInElmVelocity);
            double partial_volume = calc_determinant2x2(dxdr);
            for(int l=0; l<numOfNodeInElmVelocity; l++){
                for(int m=0; m<numOfNodeInElmPressure; m++){
                    for(int n=0; n<2; n++){
                        if(n==0){
                            PARDISO.coo_add(make_pair(element_v[ic][l], element_p[ic][m]+node.size()*2), dNdx[l][n] * Np[m]* partial_volume * gauss_weight[j] * gauss_weight[k]);
                            PARDISO.coo_add(make_pair(element_p[ic][m]+node.size()*2, element_v[ic][l]), dNdx[l][n] * Np[m]* partial_volume * gauss_weight[j] * gauss_weight[k]);
                        }
                        else if(n==1){
                            PARDISO.coo_add(make_pair(element_v[ic][l]+node.size(), element_p[ic][m]+node.size()*2), dNdx[l][n] * Np[m]* partial_volume * gauss_weight[j] * gauss_weight[k]);
                            PARDISO.coo_add(make_pair(element_p[ic][m]+node.size()*2, element_v[ic][l]+node.size()), dNdx[l][n] * Np[m]* partial_volume * gauss_weight[j] * gauss_weight[k]);
                        }
                    }
                }
            }
        }
    }
}

void STOKES_solver::main_stokes()
{
    PARDISO.initialize(node.size()*2+pressure_node.size());

    gauss_point_setting();
    //calc Kv & Kvv & Darcy
    #pragma omp parallel for
    for(int ic=0; ic<element_v.size(); ic++){
        calc_diffusive_and_Darcy_matrix(ic);
        calc_pressure_matrix(ic);
    }

    #pragma omp parallel for
    for(int ic=0;ic<node.size()*2+pressure_node.size();ic++) PARDISO.b[ic] = 0e0;

    for(int ic=0; ic<inlet_boundary_node.size(); ic++){
        for(int j=0; j<node.size()*2+pressure_node.size(); j++){
            if(PARDISO.coo_map.count(make_pair(inlet_boundary_node[ic], j))!=0){
                PARDISO.coo_insert(make_pair(inlet_boundary_node[ic], j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(inlet_boundary_node[ic]+node.size(), j))!=0){
                PARDISO.coo_insert(make_pair(inlet_boundary_node[ic]+node.size(), j),0.0);
            }
        }
        PARDISO.coo_insert(make_pair(inlet_boundary_node[ic], inlet_boundary_node[ic]), 1e0);
        PARDISO.coo_insert(make_pair(inlet_boundary_node[ic]+node.size(), inlet_boundary_node[ic]+node.size()), 1e0);
        PARDISO.b[inlet_boundary_node[ic]+node.size()] = -1e-3;
    }

    #pragma omp parallel for
    for(int j=0; j<node.size()*2+pressure_node.size(); j++){
        if(PARDISO.coo_map.count(make_pair(node.size()*2, j))!=0){
            PARDISO.coo_insert(make_pair(node.size()*2, j), 0e0);
        }
    }
    PARDISO.coo_insert(make_pair(node.size()*2, node.size()*2), 1e0);

    PARDISO.create_csr_matrix(node.size()*2+pressure_node.size());

    PARDISO.main_solve(node.size()*2+pressure_node.size(),numOfOmpThreads);

    #pragma omp parallel for
    for(int i=0; i<node.size(); i++){
        u[i] = PARDISO.x[i];
        v[i] = PARDISO.x[i+node.size()];
    }

    #pragma omp parallel for
    for(int i=0; i<pressure_node.size(); i++){
        pressure[i] = PARDISO.x[i+node.size()*2];
    }

    export_vtu_velocity("test_velocity.vtu");
    export_vtu_pressure("test_presssure.vtu");
}

void STOKES_solver::main_stokes_topology_consider_external_force(vector<vector<double>> external_force)
{
    PARDISO.initialize(node.size()*2+pressure_node.size());
    
    gauss_point_setting();
    //calc Kv & Kvv & Darcy
    #pragma omp parallel for
    for(int ic=0; ic<element_v.size(); ic++){
        calc_diffusive_and_Darcy_matrix(ic);
        calc_pressure_matrix(ic);
    }
    
    #pragma omp parallel for
    for(int ic=0;ic<node.size()*2+pressure_node.size();ic++) PARDISO.b[ic] = 0e0;


    for(int i=0;i<node.size();i++){
        for(int j=0;j<2;j++){
            PARDISO.b[i+j*node.size()] = external_force[i][j];
        }
    }
    
    for(int ic=0; ic<inlet_boundary_node.size(); ic++){
        for(int j=0; j<node.size()*2+pressure_node.size(); j++){
            if(PARDISO.coo_map.count(make_pair(inlet_boundary_node[ic], j))!=0){
                PARDISO.coo_insert(make_pair(inlet_boundary_node[ic], j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(wall_boundary_node[ic], j))!=0){
                PARDISO.coo_insert(make_pair(wall_boundary_node[ic], j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(inlet_boundary_node[ic]+node.size(), j))!=0){
                PARDISO.coo_insert(make_pair(inlet_boundary_node[ic]+node.size(), j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(wall_boundary_node[ic]+node.size(), j))!=0){
                PARDISO.coo_insert(make_pair(wall_boundary_node[ic]+node.size(), j),0.0);
            }
        }
        PARDISO.coo_insert(make_pair(inlet_boundary_node[ic], inlet_boundary_node[ic]), 1e0);
        PARDISO.coo_insert(make_pair(wall_boundary_node[ic], wall_boundary_node[ic]), 1e0);
        PARDISO.coo_insert(make_pair(inlet_boundary_node[ic]+node.size(), inlet_boundary_node[ic]+node.size()), 1e0);
        PARDISO.coo_insert(make_pair(wall_boundary_node[ic]+node.size(), wall_boundary_node[ic]+node.size()), 1e0);
        PARDISO.b[inlet_boundary_node[ic]+node.size()] = -1e-3;
    }


    for(int ic=0; ic<outlet_boundary_node.size(); ic++){
        for(int j=0; j<node.size()*2+pressure_node.size(); j++){
            if(PARDISO.coo_map.count(make_pair(outlet_boundary_node[ic]+node.size()*2, j))!=0){
                PARDISO.coo_insert(make_pair(outlet_boundary_node[ic]+node.size()*2, j),0.0);
            }
        }
        PARDISO.coo_insert(make_pair(outlet_boundary_node[ic]+node.size()*2, outlet_boundary_node[ic]+node.size()*2), 1e0);
    }



    PARDISO.coo_insert(make_pair(node.size()*2, node.size()*2), 1e0);


    PARDISO.create_csr_matrix(node.size()*2+pressure_node.size());


    PARDISO.main_solve(node.size()*2+pressure_node.size(),numOfOmpThreads);

    #pragma omp parallel for
    for(int i=0; i<node.size(); i++){
        u[i] = PARDISO.x[i];
        v[i] = PARDISO.x[i+node.size()];
    }

    #pragma omp parallel for
    for(int i=0; i<pressure_node.size(); i++){
        pressure[i] = PARDISO.x[i+node.size()*2];
    }
}

void STOKES_solver::main_stokes_topology()
{
    PARDISO.initialize(node.size()*2+pressure_node.size());
    gauss_point_setting();
    //calc Kv & Kvv & Darcy
    #pragma omp parallel for
    for(int ic=0; ic<element_v.size(); ic++){
        calc_diffusive_and_Darcy_matrix(ic);
        calc_pressure_matrix(ic);
    }

    #pragma omp parallel for
    for(int ic=0;ic<node.size()*2+pressure_node.size();ic++) PARDISO.b[ic] = 0e0;


    for(int ic=0; ic<inlet_boundary_node.size(); ic++){
        for(int j=0; j<node.size()*2+pressure_node.size(); j++){
            if(PARDISO.coo_map.count(make_pair(inlet_boundary_node[ic], j))!=0){
                PARDISO.coo_insert(make_pair(inlet_boundary_node[ic], j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(wall_boundary_node[ic], j))!=0){
                PARDISO.coo_insert(make_pair(wall_boundary_node[ic], j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(inlet_boundary_node[ic]+node.size(), j))!=0){
                PARDISO.coo_insert(make_pair(inlet_boundary_node[ic]+node.size(), j),0.0);
            }
            if(PARDISO.coo_map.count(make_pair(wall_boundary_node[ic]+node.size(), j))!=0){
                PARDISO.coo_insert(make_pair(wall_boundary_node[ic]+node.size(), j),0.0);
            }
        }
        PARDISO.coo_insert(make_pair(inlet_boundary_node[ic], inlet_boundary_node[ic]), 1e0);
        PARDISO.coo_insert(make_pair(wall_boundary_node[ic], wall_boundary_node[ic]), 1e0);
        PARDISO.coo_insert(make_pair(inlet_boundary_node[ic]+node.size(), inlet_boundary_node[ic]+node.size()), 1e0);
        PARDISO.coo_insert(make_pair(wall_boundary_node[ic]+node.size(), wall_boundary_node[ic]+node.size()), 1e0);
        PARDISO.b[inlet_boundary_node[ic]+node.size()] = -1e-3;
    }


    for(int ic=0; ic<outlet_boundary_node.size(); ic++){
        for(int j=0; j<node.size()*2+pressure_node.size(); j++){
            if(PARDISO.coo_map.count(make_pair(outlet_boundary_node[ic]+node.size()*2, j))!=0){
                PARDISO.coo_insert(make_pair(outlet_boundary_node[ic]+node.size()*2, j),0.0);
            }
        }
        PARDISO.coo_insert(make_pair(outlet_boundary_node[ic]+node.size()*2, outlet_boundary_node[ic]+node.size()*2), 1e0);
        PARDISO.b[outlet_boundary_node[ic]+node.size()*2] = 0e0;
    }


    PARDISO.coo_insert(make_pair(node.size()*2, node.size()*2), 1e0);

    PARDISO.create_csr_matrix(node.size()*2+pressure_node.size());

    PARDISO.main_solve(node.size()*2+pressure_node.size(),numOfOmpThreads);

    #pragma omp parallel for
    for(int i=0; i<node.size(); i++){
        u[i] = PARDISO.x[i];
        v[i] = PARDISO.x[i+node.size()];
    }

    #pragma omp parallel for
    for(int i=0; i<pressure_node.size(); i++){
        pressure[i] = PARDISO.x[i+node.size()*2];
    }
}