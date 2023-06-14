#include"topology.hpp"

using namespace std;

void STOKES_solver::input_domain_info_topology()
{
    string file_node = base_input_dir + "/" + "node.dat";
    string file_element = base_input_dir + "/" + "element.dat";
    string flow_path = base_input_dir + "/" + "flow_element.dat";

    int numOfElm=0;
    int numOfNode=0;
    string str;
    //element_velocity
    ifstream ifs(file_element);
    if(!ifs){
        cout << "error" << endl;
    }
    getline(ifs,str);
    numOfElm = stoi(str);
    for(int i=0; i<numOfElm; i++){
        getline(ifs,str);
        istringstream ss(str);
        vector<int> tmp_elm;
        for(int j=0; j<numOfNodeInElmVelocity; j++){
            getline(ss, str, ' ');
            tmp_elm.push_back(stoi(str));
        }
        element_v.push_back(tmp_elm);
    }
    ifs.close();
    //element_pressure
    ifs.open(file_element);
    getline(ifs,str);
    numOfElm = stoi(str);
    for(int i=0; i<numOfElm; i++){
        getline(ifs,str);
        istringstream ss(str);
        vector<int> tmp_elm;
        for(int j=0; j<numOfNodeInElmPressure; j++){
            getline(ss, str, ' ');
            tmp_elm.push_back(stoi(str));
        }
        element_p.push_back(tmp_elm);
    }
    ifs.close();
    //node_velocity
    ifs.open(file_node);
    getline(ifs,str);
    numOfNode = stoi(str);
    for(int i=0; i<numOfNode; i++){
        getline(ifs,str);
        istringstream ss(str);
        vector<double> tmp_node;
        for(int j=0; j<3; j++){
            getline(ss, str, ' ');
            tmp_node.push_back(stod(str)*1e-3);
        }
        node.push_back(tmp_node);
    }
    ifs.close();

    //phiの初期値は0.5
    phi.resize(element_v.size());
    for(int i=0; i<element_v.size(); i++){
        phi[i]=0.5;
    }
}

void STOKES_solver::input_boundary_info_topology()
{
    string inlet_wall_boundary = base_input_dir + "/" + "inlet_node.dat";
    string str;
    ifstream ifs(inlet_wall_boundary);
    while(getline(ifs,str)){
        inlet_boundary_node.push_back(stoi(str));
    }
    ifs.close();

    string outlet_wall_boundary = base_input_dir + "/" + "outlet_node.dat";
    ifs.open(outlet_wall_boundary);
    while(getline(ifs,str)){
        outlet_boundary_node.push_back(stoi(str));
    }
    ifs.close();

    string wall_boundary = base_input_dir + "/" + "wall_boundary_node.dat";
    ifs.open(wall_boundary);
    while(getline(ifs,str)){
        wall_boundary_node.push_back(stoi(str));
    }
    ifs.close();
}

double topology::calc_C()
{
    double objective_function = 0e0;
    for(int i=0; i<Stokes_main.element_v.size(); i++){
        double value=0.0;
        vector<vector<double>>  dNdr(Stokes_main.numOfNodeInElmVelocity, vector<double>(2, 0.0));
        vector<double> Np(Stokes_main.numOfNodeInElmPressure, 0.0), Nv(Stokes_main.numOfNodeInElmVelocity, 0.0);
        vector<vector<double>> dxdr(2, vector<double>(2, 0.0));
        vector<vector<double>> drdx(2, vector<double>(2, 0.0));
        vector<vector<double>> dNdx(Stokes_main.numOfNodeInElmVelocity, vector<double>(2, 0.0));
        for(int j=0; j<Stokes_main.gauss_point.size(); j++){
            for(int k=0; k<Stokes_main.gauss_point.size(); k++){
                Stokes_main.ShapeFunctionC2D8(Nv, Stokes_main.gauss_point[j], Stokes_main.gauss_point[k]);
                Stokes_main.ShapeFunctionC2D8_dNdr(dNdr, Stokes_main.gauss_point[j], Stokes_main.gauss_point[k]);
                Stokes_main.calc_dxdr(dxdr, dNdr, Stokes_main.element_v, Stokes_main.node, i, Stokes_main.numOfNodeInElmVelocity);
                Stokes_main.calc_dNdx(dNdx, drdx, dNdr, Stokes_main.numOfNodeInElmVelocity);
                double partial_volume = Stokes_main.calc_determinant2x2(dxdr);
                double dvdx[2][2] = {0e0, 0e0, 0e0, 0e0};
                for (int p = 0; p < Stokes_main.element_v[i].size(); p++){
                    dvdx[0][0] += dNdx[p][0] * Stokes_main.u[Stokes_main.element_v[i][p]];
                    dvdx[0][1] += dNdx[p][1] * Stokes_main.u[Stokes_main.element_v[i][p]];
                    dvdx[1][0] += dNdx[p][0] * Stokes_main.v[Stokes_main.element_v[i][p]];
                    dvdx[1][1] += dNdx[p][1] * Stokes_main.v[Stokes_main.element_v[i][p]];
                }

                double work=0e0;
                for(int l=0;l<2;l++){
                    for(int m=0;m<2;m++) work += Stokes_main.mu * dvdx[l][m] * (dvdx[l][m]+dvdx[m][l]) * partial_volume * Stokes_main.gauss_weight[j] * Stokes_main.gauss_weight[k];
                }

                double vel[2] = {0e0, 0e0};
                for (int l = 0; l < Stokes_main.element_v[i].size(); l++){
                    vel[0] += Nv[l] * Stokes_main.u[Stokes_main.element_v[i][l]];
                    vel[1] += Nv[l] * Stokes_main.v[Stokes_main.element_v[i][l]];
                }

                double f = Stokes_main.resistance * Stokes_main.alpha*(1e0-Stokes_main.phi[i])/(Stokes_main.alpha+Stokes_main.phi[i]);
                for (int l = 0; l < 2; l++) work += f * vel[l] * vel[l] * partial_volume * Stokes_main.gauss_weight[j] * Stokes_main.gauss_weight[k];
                objective_function += work;
            }
        }
    }
    return objective_function;
}

void topology::calc_adjoint_force()
{
    for(int i=0; i<Stokes_main.element_v.size(); i++){
        vector<vector<double>>  dNdr(Stokes_main.numOfNodeInElmVelocity, vector<double>(2, 0.0));
        vector<double> Np(Stokes_main.numOfNodeInElmPressure, 0.0), Nv(Stokes_main.numOfNodeInElmVelocity, 0.0);
        vector<vector<double>> dxdr(2, vector<double>(2, 0.0));
        vector<vector<double>> drdx(2, vector<double>(2, 0.0));
        vector<vector<double>> dNdx(Stokes_main.numOfNodeInElmVelocity, vector<double>(2, 0.0));
        for(int j=0; j<Stokes_main.gauss_point.size(); j++){
            for(int k=0; k<Stokes_main.gauss_point.size(); k++){
                Stokes_main.ShapeFunctionC2D8(Nv, Stokes_main.gauss_point[j], Stokes_main.gauss_point[k]);
                Stokes_main.ShapeFunctionC2D8_dNdr(dNdr, Stokes_main.gauss_point[j], Stokes_main.gauss_point[k]);
                Stokes_main.calc_dxdr(dxdr, dNdr, Stokes_main.element_v, Stokes_main.node, i, Stokes_main.numOfNodeInElmVelocity);
                Stokes_main.calc_dNdx(dNdx, drdx, dNdr, Stokes_main.numOfNodeInElmVelocity);
                double partial_volume = Stokes_main.calc_determinant2x2(dxdr);
                double dvdx[2][2] = {0e0, 0e0, 0e0, 0e0};
                for (int p = 0; p < Stokes_main.element_v[i].size(); p++){
                    dvdx[0][0] += dNdx[p][0] * Stokes_main.u[Stokes_main.element_v[i][p]];
                    dvdx[0][1] += dNdx[p][1] * Stokes_main.u[Stokes_main.element_v[i][p]];
                    dvdx[1][0] += dNdx[p][0] * Stokes_main.v[Stokes_main.element_v[i][p]];
                    dvdx[1][1] += dNdx[p][1] * Stokes_main.v[Stokes_main.element_v[i][p]];
                }

                for(int p=0;p<Stokes_main.element_v[i].size();p++){
                    for(int l=0; l<2; l++){
                        for(int m=0; m<2; m++){
                            adjoint_force[Stokes_main.element_v[i][p]][l] += 2e0 * Stokes_main.mu * dNdx[p][m] * (dvdx[l][m]+dvdx[m][l]) * partial_volume * Stokes_main.gauss_weight[j] * Stokes_main.gauss_weight[k];
                        }
                    }
                }

                double vel[2] = {0e0, 0e0};
                for (int l = 0; l < Stokes_main.element_v[i].size(); l++){
                    vel[0] += Nv[l] * Stokes_main.u[Stokes_main.element_v[i][l]];
                    vel[1] += Nv[l] * Stokes_main.v[Stokes_main.element_v[i][l]];
                }

                double f = Stokes_main.resistance * Stokes_main.alpha*(1e0-Stokes_main.phi[i])/(Stokes_main.alpha+Stokes_main.phi[i]);

                for(int p=0;p<Stokes_main.element_v[i].size();p++){
                    for (int l = 0; l < 2; l++){
                        adjoint_force[Stokes_main.element_v[i][p]][l] += 2e0 * f * Nv[p] * vel[l] * partial_volume * Stokes_main.gauss_weight[j] * Stokes_main.gauss_weight[k];
                    }
                }
            }
        }
    }
}

void topology::initialize()
{
    rho.resize(Stokes_main.element_v.size());
    for(int i=0; i<rho.size(); i++){
        rho[i] = 0.5;
    }
    
    adjoint_force.resize(Stokes_main.node.size());
    for(int i=0; i<adjoint_force.size(); i++){
        adjoint_force[i].resize(2);
    }
    
    dCdr.resize(Stokes_main.element_v.size());
    
    TotalVolume = 0e0;
    element_volume.resize(Stokes_main.element_v.size());

    for(int ic=0;ic<Stokes_main.element_v.size();ic++){
        element_volume[ic] = 4e0*7e-4*7e-4;
        TotalVolume += element_volume[ic];
    }

    loop = 0e0;
}

void topology::calc_sensivity()
{
    for(int i=0; i<Stokes_main.element_v.size(); i++){
        vector<vector<double>>  dNdr(Stokes_main.numOfNodeInElmVelocity, vector<double>(2, 0.0));
        vector<double> Np(Stokes_main.numOfNodeInElmPressure, 0.0), Nv(Stokes_main.numOfNodeInElmVelocity, 0.0);
        vector<vector<double>> dxdr(2, vector<double>(2, 0.0));
        vector<vector<double>> drdx(2, vector<double>(2, 0.0));
        vector<vector<double>> dNdx(Stokes_main.numOfNodeInElmVelocity, vector<double>(2, 0.0));
        for(int j=0; j<Stokes_main.gauss_point.size(); j++){
            for(int k=0; k<Stokes_main.gauss_point.size(); k++){
                Stokes_main.ShapeFunctionC2D8(Nv, Stokes_main.gauss_point[j], Stokes_main.gauss_point[k]);
                Stokes_main.ShapeFunctionC2D8_dNdr(dNdr, Stokes_main.gauss_point[j], Stokes_main.gauss_point[k]);
                Stokes_main.calc_dxdr(dxdr, dNdr, Stokes_main.element_v, Stokes_main.node, i, Stokes_main.numOfNodeInElmVelocity);
                Stokes_main.calc_dNdx(dNdx, drdx, dNdr, Stokes_main.numOfNodeInElmVelocity);
                double partial_volume = Stokes_main.calc_determinant2x2(dxdr);

                double vel[2] = {0e0, 0e0}, adj[2] = {0e0,0e0};
                for (int l = 0; l < Stokes_main.numOfNodeInElmVelocity; l++){
                    vel[0] += Nv[i] * Stokes_main.u[Stokes_main.element_v[i][l]];
                    vel[1] += Nv[i] * Stokes_main.v[Stokes_main.element_v[i][l]];

                    adj[0] += Nv[i] * Stokes_adjoint.u[Stokes_adjoint.element_v[i][l]];
                    adj[1] += Nv[i] * Stokes_adjoint.v[Stokes_adjoint.element_v[i][l]];
                }

                double f = -Stokes_main.resistance * Stokes_main.alpha * (Stokes_main.alpha+1e0) /pow(Stokes_main.alpha+0.5,2e0);

                double tmp=0e0;
                for(int l=0;l<2;l++) tmp += (vel[l]*vel[l] + vel[l]*adj[l]);
                dCdr[i] += f * tmp * partial_volume * Stokes_main.gauss_weight[j] * Stokes_main.gauss_weight[k];
            }
        }
    }
}