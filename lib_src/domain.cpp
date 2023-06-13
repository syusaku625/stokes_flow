#include"stokes.hpp"

using namespace std;

void STOKES_solver::input_domain_info()
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
    phi.resize(element_v.size());

    ifs.open(flow_path);
    while(getline(ifs, str)){
        phi[stoi(str)] = 1e0;
    }
    ifs.close();
}

void STOKES_solver::input_boundary_info()
{
    string inlet_wall_boundary = base_input_dir + "/" + "inlet_node.dat";
    string str;
    ifstream ifs(inlet_wall_boundary);
    while(getline(ifs,str)){
        inlet_boundary_node.push_back(stoi(str));
    }
    ifs.close();
}

void STOKES_solver::redifine_pressure_node_element()
{
    vector<int> memory_pressure_node_number;
    set<int> pressure_node_set;

    for(int i=0; i<element_p.size(); i++){
        for(int j=0; j<element_p[i].size(); j++){
            pressure_node_set.insert(element_p[i][j]);
        }
    }

    for(auto itr=pressure_node_set.begin(); itr!=pressure_node_set.end(); itr++){
        memory_pressure_node_number.push_back(*itr);
    }
    sort(memory_pressure_node_number.begin(), memory_pressure_node_number.end());
    pressure_node.resize(memory_pressure_node_number.size());
    
    for(int i=0; i<pressure_node.size(); i++){
        pressure_node[i].resize(2);
    }

    for(int i=0; i<memory_pressure_node_number.size(); i++){
        pressure_node[i][0]=node[memory_pressure_node_number[i]][0];
        pressure_node[i][1]=node[memory_pressure_node_number[i]][1];
        pressure_element_transform[memory_pressure_node_number[i]]=i;
    }

    for(int i=0; i<element_p.size(); i++){
        for(int j=0; j<element_p[i].size(); j++){
            element_p[i][j] = pressure_element_transform[element_p[i][j]];
        }
    }
}

void STOKES_solver::pressure_velocity_initialize()
{
    pressure.resize(pressure_node.size());
    u.resize(node.size());
    v.resize(node.size());
}