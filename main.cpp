#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<fstream>
#include<sstream>
#include<set>
#include<map>
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>
#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include<algorithm>
#include "TextParser.h"

using namespace std;

void file_read(string filename_element, string file_name_node, vector<vector<int>> &element_v, vector<vector<int>> &element_p, vector<vector<double>> &node, int numOfNodeInElmVelocity, int numOfNodeInElmPressure)
{
    int numOfElm=0;
    int numOfNode=0;
    string str;
    //element_velocity
    ifstream ifs(filename_element);
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
    ifs.open(filename_element);
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
    ifs.open(file_name_node);
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
}

void set_phi(string filename, vector<double> &phi)
{
    ifstream ifs(filename);
    string str;
    while(getline(ifs, str)){
        phi[stoi(str)] = 1e0;
    }
    ifs.close();
}

void ShapeFunctionC2D6(vector<double> &N, double g1, double g2)
{
    N[0] = g1 * (2e0 * g1 - 1e0);
    N[1] = g2 * (2e0 * g2 - 1e0);
    N[2] = (1e0 - g2 - g1) * (1e0 - 2e0 * g2 - 2e0 * g1);
    N[3] = 4e0 * g1 * g2;
    N[4] = 4e0 * g2 * (1e0 - g1- g2);
    N[5] = 4e0 * (1e0 - g1 - g2) * g1;
}

void ShapeFunctionC2D9(vector<double> &N, double g1, double g2)
{
    N[0] =  2.5e-1 * (1e+0-g1) * (1e+0-g2) * g1 * g2;
    N[1] = -2.5e-1 * (1e+0+g1) * (1e+0-g2) * g1 * g2;
    N[2] =  2.5e-1 * (1e+0+g1) * (1e+0+g2) * g1 * g2;
    N[3] = -2.5e-1 * (1e+0-g1) * (1e+0+g2) * g1 * g2;
    N[4] = -5e-1   * (1e+0-g1*g1) * (1e+0-g2) * g2;
    N[5] =  5e-1   * (1e+0+g1) * (1e+0-g2*g2) * g1;
    N[6] =  5e-1   * (1e+0-g1*g1) * (1e+0+g2) * g2;
    N[7] = -5e-1   * (1e+0-g1) * (1e+0-g2*g2) * g1;
    N[8] =           (1e+0-g1*g1) * (1e+0-g2*g2);
}

void ShapeFunctionC2D3(vector<double> &N, double g1, double g2)
{
    N[0] = g1;
    N[1] = g2;
    N[2] = 1e0 - g1 - g2;
}

void ShapeFunctionC2D4(vector<double> &N, double g1, double g2)
{
    N[0] = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
    N[1] = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
    N[2] = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
    N[3] = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
}

void ShapeFunctionC2D9_dNdr(vector<vector<double>> &dNdr, double g1, double g2)
{
    dNdr[0][0] = 2.5e-1  * (1e+0-g2) * g2 * (1e+0 - 2e+0*g1);
    dNdr[0][1] = 2.5e-1  * (1e+0-g1) * g1 * (1e+0 - 2e+0*g2);
    dNdr[1][0] = -2.5e-1 * (1e+0-g2) * g2 * (1e+0 + 2e+0*g1);
    dNdr[1][1] = -2.5e-1 * (1e+0+g1) * g1 * (1e+0 - 2e+0*g2);
    dNdr[2][0] = 2.5e-1  * (1e+0+g2) * g2 * (1e+0 + 2e+0*g1);
    dNdr[2][1] = 2.5e-1  * (1e+0+g1) * g1 * (1e+0 + 2e+0*g2);
    dNdr[3][0] = -2.5e-1 * (1e+0+g2) * g2 * (1e+0 - 2e+0*g1);
    dNdr[3][1] = -2.5e-1 * (1e+0-g1) * g1 * (1e+0 + 2e+0*g2);
    dNdr[4][0] = g1 * (1e+0-g2) * g2;
    dNdr[4][1] = -5e-1 * (1e+0-g1*g1) * (1e+0-2e+0*g2);
    dNdr[5][0] = 5e-1 * (1e+0-g2*g2) * (1e+0+2e+0*g1);
    dNdr[5][1] = -g2 * (1e+0+g1) * g1;
    dNdr[6][0] = -g1 * (1e+0+g2) * g2;
    dNdr[6][1] = 5e-1    * (1e+0-g1*g1) * (1e+0+2e+0*g2);
    dNdr[7][0] = -5e-1   * (1e+0-g2*g2) * (1e+0-2e+0*g1);
    dNdr[7][1] = g2 * (1e+0-g1) * g1;
    dNdr[8][0] = -2e+0 * g1 * (1e+0-g2*g2);
    dNdr[8][1] = -2e+0 * g2 * (1e+0-g1*g1);
}

void ShapeFunctionC2D6_dNdr(vector<vector<double>> &dNdr, double g1, double g2)
{
    dNdr[0][0] = 4e0*g1-1e0;
    dNdr[0][1] = 0e0;
    dNdr[1][0] = 0e0;
    dNdr[1][1] = 4e0*g2-1e0;
    dNdr[2][0] = -3e0 + 4e0 * g1 + 4e0 * g2;
    dNdr[2][1] = -3e0 + 4e0 * g1 + 4e0 * g2;
    dNdr[3][0] = 4e0*g2;
    dNdr[3][1] = 4e0*g1;
    dNdr[4][0] = -4e0*g2;
    dNdr[4][1] = 4e0 * (1e0 - g1 - 2e0 * g2);
    dNdr[5][0] = 4e0 * (1e0 - g2 - 2e0 * g1);
    dNdr[5][1] = -4e0*g1;
}

void ShapeFunctionC2D3_dNdr(vector<vector<double>> &dNdr, double g1, double g2)
{
    dNdr[0][0] = 1e0;
    dNdr[0][1] = 0e0;
    dNdr[1][0] = 0e0;
    dNdr[1][1] = 1e0;
    dNdr[2][0] = -1e0;
    dNdr[2][1] = -1e0;
}

void calc_dxdr(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr, vector<vector<int>> element, vector<vector<double>> x, int ic, int numOfNodeInElm)
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

void calc_dxdr_C2D3(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr, vector<vector<int>> element, vector<vector<double>> x, int ic)
{
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            dxdr[i][j] = 0e0;
            for (int k = 0; k < 3; k++){
                dxdr[i][j] += dNdr[k][j] * x[element[ic][k]][i];
            }
        }
    }
}

void calc_dNdx(vector<vector<double>> &dNdx, vector<vector<double>> drdx, vector<vector<double>> &dNdr, int numOfNodeInElm)
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

void calc_dNdx_C2D3(vector<vector<double>> &dNdx, vector<vector<double>> drdx, vector<vector<double>> &dNdr)
{
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 2; j++){
            dNdx[i][j] = 0e0;
            for (int k = 0; k < 2; k++){
                dNdx[i][j] += dNdr[i][k] * drdx[k][j];
            }
        }
    }
}

double calc_determinant2x2(std::vector<std::vector<double>> A)
{
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];
}

void calc_inverse_matrix2x2(std::vector<std::vector<double>> &drdx, std::vector<std::vector<double>> dxdr)
{
    double inverse=calc_determinant2x2(dxdr);
    drdx[0][0]=1e0/inverse*dxdr[1][1];
    drdx[1][1]=1e0/inverse*dxdr[0][0];
    drdx[0][1]=1e0/inverse*(dxdr[0][1])*(-1);
    drdx[1][0]=1e0/inverse*(dxdr[1][0])*(-1);
}

void redefine_pressure_node_element(vector<vector<double>> node, vector<vector<int>> &element_p, vector<vector<double>> &pressure_node, map<int, int> &pressure_element_transform)
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

void export_vtu_velocity(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> u, vector<double> v)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", int(node.size()), int(element.size()));
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < static_cast<int>(element[i].size()); j++) fprintf(fp, "%d ", int(element[i][j]));
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += static_cast<int>(element[i].size());
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    if(element[i].size()==6) fprintf(fp, "%d\n", 22);
    if(element[i].size()==9) fprintf(fp, "%d\n", 23);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity[m/s]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = u[ic];
      data_d[num+1]   = v[ic];
      data_d[num+2]   = 0e0;
      num+=3;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  
  delete[] data_d;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void export_vtu_pressure(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> p)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", int(node.size()), int(element.size()));
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < static_cast<int>(element[i].size()); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += static_cast<int>(element[i].size());
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    if(element[i].size()==3) fprintf(fp, "%d\n", 5);
    if(element[i].size()==4) fprintf(fp, "%d\n", 9);
  }
    
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size();
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = p[ic];
      num++;
  }
  size=sizeof(double)*node.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  
  delete[] data_d;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

double  return_max_edge_length(vector<vector<double>> node, vector<vector<int>> element, int ic)
{
    double L1 = pow(pow((node[element[ic][0]][0] - node[element[ic][1]][0]),2.0) + pow((node[element[ic][0]][1] - node[element[ic][1]][1]),2.0), 0.5);
    double L2 = pow(pow((node[element[ic][1]][0] - node[element[ic][1]][2]),2.0) + pow((node[element[ic][1]][1] - node[element[ic][2]][1]),2.0), 0.5); 
    double L3 = pow(pow((node[element[ic][2]][0] - node[element[ic][1]][0]),2.0) + pow((node[element[ic][2]][1] - node[element[ic][0]][1]),2.0), 0.5);
    return max(max(L1, L2), L3);
}

typedef Eigen::Triplet<double> T;

int main(int argc,char *argv[])
{
    TextParser tp;
    std::string input_file = argv[1];
    int ierror;
    if ((ierror = tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }
    string base_label = "/Domain";
    string label = base_label + "/base_folder";
    string base_input_dir;
    double alpha, resistance, mu;
    vector<vector<int>> element_v, element_p;
    vector<vector<double>> node, pressure_node;
    map<int, int> pressure_element_transform;
    vector<double> phi;
    string geometry;
    vector<double> pressure;
    vector<double> u, v;

    if ( !tp.getInspectedValue(label,base_input_dir)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    int numofNodeinElm_velocity, numofNodeinElm_pressure;
    label = base_label + "/numOfNodeInElmVelocity";
    if ( !tp.getInspectedValue(label,numofNodeinElm_velocity)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    label = base_label + "/numOfNodeInElmPressure";
    if ( !tp.getInspectedValue(label,numofNodeinElm_pressure)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    label = base_label + "/geometry";
    if ( !tp.getInspectedValue(label,geometry)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    base_label = "/Parameter";
    label = base_label + "/mu";
    if ( !tp.getInspectedValue(label,mu)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    if(geometry == "quad"){
        label = base_label + "/alpha";
        if ( !tp.getInspectedValue(label,alpha)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/resistance";
        if ( !tp.getInspectedValue(label,resistance)){
            cout << label << " is not set" << endl;
            exit(0);
        }
    }

    string file_node = base_input_dir + "/" + "node.dat";
    string file_element = base_input_dir + "/" + "element.dat";
    string flow_path = base_input_dir + "/" + "flow_element.dat";

    file_read(file_element, file_node, element_v, element_p, node, numofNodeinElm_velocity, numofNodeinElm_pressure);
    
    if(geometry=="quad"){
        phi.resize(element_v.size());
        set_phi(flow_path, phi);
    }

    redefine_pressure_node_element(node, element_p, pressure_node, pressure_element_transform);
    pressure.resize(pressure_node.size());
    u.resize(node.size());
    v.resize(node.size());

    vector<double> gauss_point, gauss_weight;
    if (geometry == "triangle") {
        gauss_point.resize(3); gauss_weight.resize(3);
        gauss_point[0] = 2e0/3e0; gauss_point[1] = 1e0/6e0; gauss_point[2] = 1e0/6e0;
        gauss_weight[0] = 1e0/3e0; gauss_weight[1] = 1e0/3e0; gauss_weight[2] = 1e0/3e0;
    }
    else if (geometry == "quad"){
        gauss_point.resize(3); gauss_weight.resize(3);
        gauss_point[0] = -0.774596669241483; gauss_point[1] = 0e0; gauss_point[2] = 0.774596669241483;
        gauss_weight[0] = 0.555555555555555; gauss_weight[1] = 0.888888888888888; gauss_weight[2] = 0.555555555555555;
    }
    
    vector<vector<double>> Kv(node.size(), vector<double>(node.size(), 0.0));
    vector<vector<vector<double>>> Kp(node.size(), vector<vector<double>>(pressure_node.size(), vector<double>(2, 0.0)));
    vector<vector<double>> Darcy(node.size(), vector<double>(node.size(), 0.0));
    vector<vector<double>> global_matrix(node.size()*2+pressure_node.size(), vector<double>(node.size()*2+pressure_node.size(), 0.0));
    vector<double> b(node.size()*2+pressure_node.size(), 0.0);

    //calc Kv & Kvv & Darcy
    for(int i=0; i<element_v.size(); i++){
        vector<vector<double>> dxdr(2, vector<double>(2, 0.0)), drdx(2, vector<double>(2, 0.0)), dNdr(numofNodeinElm_velocity, vector<double>(2, 0.0)), dNdx(numofNodeinElm_velocity, vector<double>(2, 0.0));
        vector<double> Np(numofNodeinElm_pressure, 0.0), Nv(numofNodeinElm_velocity, 0.0);
        double sum_volume = 0.0;
        for(int j=0; j<gauss_point.size(); j++){
            for(int k=0; k<gauss_point.size(); k++){
                if (geometry == "triangle"){
                    ShapeFunctionC2D6_dNdr(dNdr, gauss_point[j], gauss_point[k]);
                    ShapeFunctionC2D6(Nv, gauss_point[j], gauss_point[k]);
                    ShapeFunctionC2D3(Np, gauss_point[j], gauss_point[k]);
                    calc_dxdr(dxdr, dNdr, element_v, node, i, numofNodeinElm_velocity);
                    calc_inverse_matrix2x2(drdx,dxdr);
                    calc_dNdx(dNdx, drdx, dNdr, numofNodeinElm_velocity);
                }
                else if (geometry == "quad"){
                    ShapeFunctionC2D9_dNdr(dNdr, gauss_point[j], gauss_point[k]);
                    ShapeFunctionC2D9(Nv, gauss_point[j], gauss_point[k]);
                    ShapeFunctionC2D4(Np, gauss_point[j], gauss_point[k]);
                    calc_dxdr(dxdr, dNdr, element_v, node, i, numofNodeinElm_velocity);
                    calc_inverse_matrix2x2(drdx,dxdr);
                    calc_dNdx(dNdx, drdx, dNdr, numofNodeinElm_velocity);
                }
                double partial_volume = calc_determinant2x2(dxdr);
                sum_volume += partial_volume * gauss_weight[j] * gauss_weight[k];
                //6x6 matrix
                for(int l=0; l<numofNodeinElm_velocity; l++){
                    for(int m=0; m<numofNodeinElm_velocity; m++){
                        for(int n=0; n<2; n++){
                            Kv[element_v[i][l]][element_v[i][m]] -= mu * dNdx[l][n]*dNdx[m][n] * partial_volume * gauss_weight[j] * gauss_weight[k];
                        }
                    }
                }
                //6x3 matrix
                for(int l=0; l<numofNodeinElm_velocity; l++){
                    for(int m=0; m<numofNodeinElm_pressure; m++){
                        for(int n=0; n<2; n++){
                            Kp[element_v[i][l]][element_p[i][m]][n] += dNdx[l][n] * Np[m]* partial_volume * gauss_weight[j] * gauss_weight[k];
                        }
                    }
                }
                if(geometry == "quad"){
                    //Darcy matrix
                    for(int l=0; l<numofNodeinElm_velocity; l++){
                        for(int m=0; m<numofNodeinElm_velocity; m++){
                            Darcy[element_v[i][l]][element_v[i][m]] -= resistance * alpha * (1e0 - phi[i]) / (alpha + phi[i]) * Nv[l] * Nv[m] * partial_volume * gauss_weight[j] * gauss_weight[k]; 
                        }
                    }
                }
            }
        }
    }

    //set Kv& Darcy
    for(int i=0; i<node.size(); i++){
        for(int j=0; j<node.size(); j++){
            global_matrix[i][j] = Kv[i][j];
            global_matrix[i+node.size()][j+node.size()] = Kv[i][j];
            if(geometry=="quad"){
                global_matrix[i][j] += Darcy[i][j];
                global_matrix[i+node.size()][j+node.size()] += Darcy[i][j];
            }
        }
    }
    //set Kp
    for(int i=0; i<node.size(); i++){
        for(int j=0; j<pressure_node.size(); j++){
            global_matrix[i][j+2*node.size()] = Kp[i][j][0];
            global_matrix[i+node.size()][j+2*node.size()] = Kp[i][j][1];
        }
    }
    //set Kvv u
    for(int i=0; i<pressure_node.size(); i++){
        for(int j=0; j<node.size(); j++){
            global_matrix[i+2*node.size()][j] = Kp[j][i][0];
            global_matrix[i+2*node.size()][j+node.size()] = Kp[j][i][1];
        }
    }

    string str;
    vector<int> wall_boundary_node, inlet_boundary_node, outlet_boundary_node;
    string inlet_wall_boundary = base_input_dir + "/" + "inlet_node.dat";
    ifstream ifs(inlet_wall_boundary);
    while(getline(ifs,str)){
        inlet_boundary_node.push_back(stoi(str));
    }
    ifs.close();
    if(geometry=="quad"){
        for(int i=0; i<inlet_boundary_node.size(); i++){
            for(int j=0; j<global_matrix.size(); j++){
                global_matrix[inlet_boundary_node[i]][j] = 0.0; //u
                global_matrix[inlet_boundary_node[i]+node.size()][j] = 0.0; //v
            }
            global_matrix[inlet_boundary_node[i]][inlet_boundary_node[i]] = 1.0;
            global_matrix[inlet_boundary_node[i]+node.size()][inlet_boundary_node[i]+node.size()] = 1.0;
            b[inlet_boundary_node[i]+node.size()] = -1e-3;
        }

        for(int j=0; j<global_matrix.size(); j++){
            global_matrix[node.size()*2][j] = 0.0; //p
        }
        global_matrix[node.size()*2][node.size()*2] = 1.0;
    }

    if(geometry=="triangle"){
        string outlet_wall_boundary = base_input_dir + "/" + "pressure_wall_node.dat";
        string wall_boundary = base_input_dir + "/" + "wall_node.dat";

        ifstream ifs(wall_boundary);
        while(getline(ifs,str)){
            wall_boundary_node.push_back(stoi(str));
        }
        ifs.close();
        ifs.open(outlet_wall_boundary);
        while(getline(ifs,str)){
            outlet_boundary_node.push_back(stoi(str));
        }
        ifs.close();
        for(int i=0; i<wall_boundary_node.size(); i++){
            for(int j=0; j<global_matrix.size(); j++){
                global_matrix[wall_boundary_node[i]][j] = 0.0; //u
                global_matrix[wall_boundary_node[i]+node.size()][j] = 0.0; //v
            }
            global_matrix[wall_boundary_node[i]][wall_boundary_node[i]] = 1.0;
            global_matrix[wall_boundary_node[i]+node.size()][wall_boundary_node[i]+node.size()] = 1.0;
        }

        for(int i=0; i<inlet_boundary_node.size(); i++){
            for(int j=0; j<global_matrix.size(); j++){
                global_matrix[inlet_boundary_node[i]][j] = 0.0; //u
                global_matrix[inlet_boundary_node[i]+node.size()][j] = 0.0; //v
            }
            global_matrix[inlet_boundary_node[i]][inlet_boundary_node[i]] = 1.0;
            global_matrix[inlet_boundary_node[i]+node.size()][inlet_boundary_node[i]+node.size()] = 1.0;
            b[inlet_boundary_node[i]] = 1e-3;
        }
        for(int i=0; i<outlet_boundary_node.size(); i++){
            for(int j=0; j<global_matrix.size(); j++){
                global_matrix[outlet_boundary_node[i]+node.size()*2][j] = 0.0; //p
            }
            global_matrix[outlet_boundary_node[i]+node.size()*2][outlet_boundary_node[i]+node.size()*2] = 1.0;
        }
    }


    std::vector<T> tripletVec;
    for(int i=0; i<global_matrix.size(); i++){
        for(int j=0; j<global_matrix[i].size(); j++){
            if(fabs(global_matrix[i][j])>1e-12){
                tripletVec.push_back(T(i,j,global_matrix[i][j]));
            }
        }
    }
    
    Eigen::SparseMatrix<double> M(global_matrix.size(),global_matrix.size());
    M.setFromTriplets(tripletVec.begin(), tripletVec.end());
    Eigen::SparseMatrix<double> b_eigen(b.size(),1);
    std::vector<T> tripletList;
    for(int i=0; i<b.size(); i++){
        if(fabs(b[i])>1e-12){
            tripletList.push_back(T(i,0,b[i]));
        }
    }
    b_eigen.setFromTriplets(tripletList.begin(), tripletList.end());
    Eigen::VectorXd x;  // 解のベクトル
    Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver2;  // solverオブジェクトを構築する。

    solver2.compute(M);
    if( solver2.info() != Eigen::Success ) {
        std::cerr << "decomposition failed" << std::endl;
        exit(1);
    }
    x = solver2.solve(b_eigen);
    if( solver2.info() != Eigen::Success ) {
        std::cerr << "solving failed" << std::endl;
        exit(1);
    }
    
    for(int i=0; i<node.size(); i++){
        u[i] = x(i);
        v[i] = x(i+node.size());
    }
    for(int i=0; i<pressure_node.size(); i++){
        pressure[i] = x(i+node.size()*2);
    }

    export_vtu_velocity("test_velocity.vtu", element_v, node, u, v);
    export_vtu_pressure("test_presssure.vtu", element_p, pressure_node, pressure);
}