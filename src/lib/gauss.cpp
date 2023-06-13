#include"stokes.hpp"

using namespace std;

void STOKES_solver::gauss_point_setting()
{
    gauss_point.resize(3); gauss_weight.resize(3);
    gauss_point[0] = -0.774596669241483; gauss_point[1] = 0e0; gauss_point[2] = 0.774596669241483;
    gauss_weight[0] = 0.555555555555555; gauss_weight[1] = 0.888888888888888; gauss_weight[2] = 0.555555555555555;   
}