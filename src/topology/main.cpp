#include"topology.hpp"
#include"MMA.hpp"

using namespace std;

int main(int argc,char *argv[])
{
    topology T_S;
    MMA mma;
    
    std::string input_file = argv[1];
    int ierror;

    if ((ierror = T_S.tp.read(input_file)) != TP_NO_ERROR) {
        printf("\tError at reading '%s' file\n", input_file.c_str());
        return 1;
    }

    string base_label = "/Domain";
    string label = base_label + "/loop_max";

    if ((ierror = T_S.Stokes_main.tp.read(input_file)) != TP_NO_ERROR) {
        printf("\tError at reading '%s' file\n", input_file.c_str());
        return 1;
    }
    if ((ierror = T_S.Stokes_adjoint.tp.read(input_file)) != TP_NO_ERROR) {
        printf("\tError at reading '%s' file\n", input_file.c_str());
        return 1;
    }

    T_S.Stokes_main.input_parameter_from_tp_file();
    T_S.Stokes_main.input_domain_info_topology();
    T_S.Stokes_main.input_boundary_info_topology();
    T_S.Stokes_main.redifine_pressure_node_element();
    T_S.Stokes_main.pressure_velocity_initialize();

    T_S.Stokes_adjoint.input_parameter_from_tp_file();
    T_S.Stokes_adjoint.input_domain_info_topology();
    T_S.Stokes_adjoint.input_boundary_info_topology();
    T_S.Stokes_adjoint.redifine_pressure_node_element();
    T_S.Stokes_adjoint.pressure_velocity_initialize();

    T_S.initialize();

    mma.simp_usingMMA(T_S);
}