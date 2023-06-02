#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<fstream>
#include<sstream>
#include<set>
#include<map>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include<algorithm>

using namespace std;

void file_read(string filename_element, string file_name_node, vector<vector<int>> &element_v, vector<vector<int>> &element_p, vector<vector<double>> &node)
{
    int numOfElm=0;
    int numOfNode=0;
    string str;
    //element_velocity
    ifstream ifs(filename_element);
    getline(ifs,str);
    numOfElm = stoi(str);
    for(int i=0; i<numOfElm; i++){
        getline(ifs,str);
        istringstream ss(str);
        vector<int> tmp_elm;
        for(int j=0; j<6; j++){
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
        for(int j=0; j<3; j++){
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

void ShapeFunctionC2D6(vector<double> &N, double g1, double g2)
{
    N[0] = g1 * (2e0 * g1 - 1e0);
    N[1] = g2 * (2e0 * g2 - 1e0);
    N[2] = (1e0 - g2 - g1) * (1e0 - 2e0 * g2 - 2e0 * g1);
    N[3] = 4e0 * g1 * g2;
    N[4] = 4e0 * g2 * (1e0 - g1- g2);
    N[5] = 4e0 * (1e0 - g1 - g2) * g1;
}

void ShapeFunctionC2D3(vector<double> &N, double g1, double g2)
{
    N[0] = g1;
    N[1] = g2;
    N[2] = 1e0 - g1 - g2;
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

void calc_dxdr(std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr, vector<vector<int>> element, vector<vector<double>> x, int ic)
{
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            dxdr[i][j] = 0e0;
            for (int k = 0; k < 6; k++){
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

void calc_dNdx(vector<vector<double>> &dNdx, vector<vector<double>> drdx, vector<vector<double>> &dNdr)
{
    for (int i = 0; i < 6; i++){
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
    drdx[0][0]=1e0/fabs(inverse)*dxdr[1][1];
    drdx[1][1]=1e0/fabs(inverse)*dxdr[0][0];
    drdx[0][1]=1e0/fabs(inverse)*(dxdr[0][1])*(-1);
    drdx[1][0]=1e0/fabs(inverse)*(dxdr[1][0])*(-1);
}

void Jacobi_method(std::vector<double> &p, vector<double> &R, vector<vector<double>> &G, double convergence)
{
    vector<double> tmp_p;
    int NN=p.size();
    while(1){
        tmp_p.resize(NN);
        for (int i = 0; i < NN; i++)
        {
            tmp_p[i]=R[i];
            for (int j = 0; j < NN; j++){
                if(i!=j){
                    tmp_p[i] -= G[i][j] * tmp_p[j];
                }
            }
            tmp_p[i] /= G[i][i];
            //cout << G[i][i] << endl;
            //cout << tmp_p[i] << endl;
        }
        exit(1);
        double err = 0.0;
        for (int i = 0; i < NN; i++){
            err += fabs(tmp_p[i] - p[i]);
            p[i] = tmp_p[i];
        }
        //cout << err << endl;
        if (err < convergence){
            break;
        }
    }
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

void export_vtu_quadratic_triangle(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> u, vector<double> v)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
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
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 22);
    
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

void export_vtu_triangle(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> p)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
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
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 5);
    
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

int main()
{
    string file_node = "test_quadratic_node.dat";
    string file_element = "test_quadratic_element.dat";
    vector<vector<int>> element_v, element_p;
    vector<vector<double>> node;
    file_read(file_element, file_node, element_v, element_p, node);

    vector<vector<double>> pressure_node;
    map<int, int> pressure_element_transform;
    redefine_pressure_node_element(node, element_p, pressure_node, pressure_element_transform);

    vector<double> pressure(pressure_node.size());
    vector<double> u(node.size()), v(node.size());

    vector<double> gauss_point = {-0.861135311594053, -0.339981043584856, 0.339981043584856, 0.861135311594053};
    vector<double> gauss_weight = {0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};
    
    //calc Kv & Kvv
    vector<vector<double>> Kv(node.size(), vector<double>(node.size(), 0.0));
    vector<vector<vector<double>>> Kp(node.size(), vector<vector<double>>(pressure_node.size(), vector<double>(2, 0.0)));
    for(int i=0; i<element_v.size(); i++){
        vector<vector<double>> dxdr(2, vector<double>(2)), drdx(2, vector<double>(2)), dNdr(6, vector<double>(2)), dNdx(6, vector<double>(2));
        vector<double> N(3);
        double sum_volume = 0.0;
        for(int j=0; j<gauss_point.size(); j++){
            for(int k=0; k<gauss_point.size(); k++){
                ShapeFunctionC2D6_dNdr(dNdr, gauss_point[j], gauss_point[k]);
                calc_dxdr(dxdr, dNdr, element_v, node, i);
                calc_inverse_matrix2x2(drdx,dxdr);
                calc_dNdx(dNdx, drdx, dNdr);
                ShapeFunctionC2D3(N, gauss_point[j], gauss_point[k]);
                double partial_volume = calc_determinant2x2(dxdr);
                sum_volume += partial_volume * gauss_weight[j] * gauss_weight[k];
                //6x6 matrix
                for(int l=0; l<6; l++){
                    for(int m=0; m<6; m++){
                        for(int n=0; n<2; n++){
                            Kv[element_v[i][l]][element_v[i][m]] -= 1e-3 * dNdx[l][n]*dNdx[m][n] * partial_volume * gauss_weight[j] * gauss_weight[k];
                        }
                    }
                }
                //6x3 matrix
                for(int l=0; l<6; l++){
                    for(int m=0; m<3; m++){
                        for(int n=0; n<2; n++){
                            Kp[element_v[i][l]][element_p[i][m]][n] += dNdx[l][n] * N[m]* partial_volume * gauss_weight[j] * gauss_weight[k];
                        }
                    }
                }
            }
        }
        //cout << sum_volume << endl;
    }

    //calc PSPG term
    //vector<vector<double>> Kpspg(pressure_node.size(), vector<double>(pressure_node.size(), 0.0));
    //for(int i=0; i<element_p.size(); i++){
    //    vector<vector<double>> dxdr(2, vector<double>(2)), drdx(2, vector<double>(2)), dNdr(3, vector<double>(2)), dNdx(3, vector<double>(2));
    //    vector<double> N(3);
    //    for(int j=0; j<gauss_point.size(); j++){
    //        for(int k=0; k<gauss_point.size(); k++){
    //            ShapeFunctionC2D3_dNdr(dNdr, gauss_point[j], gauss_point[k]);
    //            calc_dxdr_C2D3(dxdr, dNdr, element_p, pressure_node, i);
    //            calc_inverse_matrix2x2(drdx,dxdr);
    //            calc_dNdx_C2D3(dNdx, drdx, dNdr);
    //            ShapeFunctionC2D3(N, gauss_point[j], gauss_point[k]);
    //            double partial_volume = calc_determinant2x2(dxdr);
    //            double he = return_max_edge_length(pressure_node, element_p, i);
    //            //cout << he << endl;
    //            double tau = pow(pow((2e0*1e-6/he),2.0) + (4e0/pow(he,2.0)),-0.5);
    //            //3x3 matrix
    //            for(int l=0; l<3; l++){
    //                for(int m=0; m<3; m++){
    //                    for(int n=0; n<2; n++){
    //                        Kpspg[element_p[i][l]][element_p[i][m]] += tau * (dNdx[l][n] * dNdx[m][n]); 
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

    vector<vector<double>> global_matrix(node.size()*2+pressure_node.size(), vector<double>(node.size()*2+pressure_node.size(), 0.0));
    vector<double> b(node.size()*2+pressure_node.size(), 0.0);
    //set Kv
    for(int i=0; i<node.size(); i++){
        for(int j=0; j<node.size(); j++){
            global_matrix[i][j] = Kv[i][j];
            global_matrix[i+node.size()][j+node.size()] = Kv[i][j];
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
    //set Kpspg
    //for(int i=0; i<pressure_node.size(); i++){
    //    for(int j=0; j<pressure_node.size(); j++){
    //        global_matrix[i+2*node.size()][j+2*node.size()] = Kpspg[i][j];
    //    }
    //}

    ofstream test("test.csv");
    for(int i=0; i<global_matrix.size(); i++){
        for(int j=0; j<global_matrix[i].size(); j++){
            test << global_matrix[i][j] << ",";
        }
        test << endl;
    }
    string str;
    vector<int> wall_boundary_node, inlet_boundary_node, outlet_boundary_node;
    ifstream ifs("upper_wall_node.dat");
    while(getline(ifs,str)){
        wall_boundary_node.push_back(stoi(str));
    }
    ifs.close();
    ifs.open("lower_wall_node.dat");
    while(getline(ifs,str)){
        wall_boundary_node.push_back(stoi(str));
    }
    ifs.close();
    ifs.open("inlet_wall_boundary.dat");
    while(getline(ifs,str)){
        inlet_boundary_node.push_back(stoi(str));
    }
    ifs.close();
    ifs.open("pressure_wall_node.dat");
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
    std::vector<T> tripletVecb;
    for(int i=0; i<b.size(); i++){
        if(fabs(b[i])>1e-12){
            tripletVecb.push_back(T(i,0,b[i]));
        }
    }
    b_eigen.setFromTriplets(tripletVecb.begin(), tripletVecb.end());
    Eigen::VectorXd x;  // 解のベクトル
    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;  // solverオブジェクトを構築する。
    solver.compute(M);
    if( solver.info() != Eigen::Success ) {
        std::cerr << "decomposition failed" << std::endl;
    }
    x = solver.solve(b_eigen);
    if( solver.info() != Eigen::Success ) {
        std::cerr << "solving failed" << std::endl;
    }
    ofstream ofs("ans.dat");
    for(int i=0; i<node.size(); i++){
        u[i] = x(i);
        v[i] = x(i+node.size());
    }
    for(int i=0; i<pressure_node.size(); i++){
        pressure[i] = x(i+node.size()*2);
    }

    export_vtu_quadratic_triangle("test_velocity.vtu", element_v, node, u, v);
    export_vtu_triangle("test_presssure.vtu", element_p, pressure_node, pressure);
}