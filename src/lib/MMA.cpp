#include "MMA.hpp"

using namespace std;
double MMA::return_objective_function(const vector<double> &rho_MMA, vector<double> &grad, void *my_func_data)
{
    topology *topo_ptr  = reinterpret_cast<topology*>(my_func_data);

    for(int i=0; i<topo_ptr->Stokes_main.element_v.size(); i++){
        topo_ptr->dCdr[i] = 0e0;
    }

    for(int i=0; i<topo_ptr->Stokes_main.element_v.size(); i++){
        topo_ptr->Stokes_main.phi[i] = rho_MMA[i];
        topo_ptr->Stokes_adjoint.phi[i] = rho_MMA[i];
    }
    
    topo_ptr->Stokes_main.main_stokes_topology();
    
    //calc C
    double objective_function = topo_ptr->calc_C();

    //calc adjoint force
    topo_ptr->calc_adjoint_force();

    topo_ptr->Stokes_adjoint.main_stokes_topology_consider_external_force(topo_ptr->adjoint_force);

    //calc dLdx
    topo_ptr->calc_sensivity();


    #pragma omp parallel for
    for(int i=0; i<topo_ptr->Stokes_main.element_v.size(); i++){  //initalize
        grad[i] = 0e0;
    }

    //no filtered grad
    for(int e=0;e<topo_ptr->Stokes_main.element_v.size();e++){
        grad[e] = topo_ptr->dCdr[e];
    }

    topo_ptr->loop++;
    string output = "simp/simp_" + to_string(topo_ptr->loop)+ ".vtu";
    cout << "loop = " << topo_ptr->loop << endl;
    topo_ptr->Stokes_main.export_vtu_velocity(output);
    return objective_function;
}

void MMA::simp_usingMMA(topology &simp)
{
    simp.loop = 0;
    nlopt::opt opt(nlopt::LD_MMA, simp.Stokes_main.element_v.size());
    
    opt.set_lower_bounds(1e-5);
    opt.set_upper_bounds(1e0);

    opt.set_min_objective(MMA::return_objective_function, &simp);

    opt.add_inequality_constraint(MMA::constraint_volume, &simp, 1e-8);
    
    opt.set_ftol_rel(1e-6);
    opt.set_maxeval(500); //for erro code: n回の反復で最適化を終了

    double minimize_objfunction;
    
    try{
        nlopt::result result = opt.optimize(simp.rho, minimize_objfunction);
        cout << "---------optimization loop completed!----------" << endl;
        cout << "optimal_objfunc=" << minimize_objfunction << endl;
    }
    catch(exception &e) {
        cout << "nlopt failed: " << e.what() << endl;
    }
}

double MMA::constraint_volume(const vector<double> &rho_MMA, vector<double> &grad, void *my_func_data)
{
    topology *topo_ptr  = reinterpret_cast<topology*>(my_func_data);
    topo_ptr->currentVolume = 0e0;

    for(int i=0;i<topo_ptr->Stokes_main.element_v.size();i++){
        topo_ptr->currentVolume+=rho_MMA[i]*topo_ptr->element_volume[i];
    }

  //cout << " constfunc = "<< left << setw(7) << simp_ptr->currentVolume << endl;

  if (!grad.empty()) {
    #pragma omp parallel for
    for(int i=0; i<topo_ptr->Stokes_main.element_v.size(); i++){  //initalize
        grad[i] = 0e0;
    }

    //no filtered grad
    for(int e=0;e<topo_ptr->Stokes_main.element_v.size();e++){
        grad[e] = topo_ptr->element_volume[e];
    }
  }


  //FILE *fp;
  //fp=fopen("volume.log","a");
  //fprintf(fp,"%d %e %e\n",simp_ptr->SIMP_loop,simp_ptr->currentVolume, simp_ptr->currentVolume / simp_ptr->TotalVolume);
  //fclose(fp);
  double Volume_max = topo_ptr->TotalVolume*0.3;
  return (topo_ptr->currentVolume - Volume_max);
}