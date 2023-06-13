#include"stokes.hpp"

using namespace std;

int main(int argc,char *argv[])
{
    STOKES_solver Stokes;
    std::string input_file = argv[1];
    int ierror;
    if ((ierror = Stokes.tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }

    Stokes.input_parameter_from_tp_file();
    Stokes.input_domain_info();
    Stokes.input_boundary_info();

    Stokes.redifine_pressure_node_element();

    Stokes.pressure_velocity_initialize();

    Stokes.main_stokes();
}