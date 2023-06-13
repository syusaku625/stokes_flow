#include"stokes.hpp"

using namespace std;

void STOKES_solver::input_parameter_from_tp_file()
{
    string base_label = "/Domain";
    string label = base_label + "/base_folder";
    
    if ( !tp.getInspectedValue(label,base_input_dir)){
        cout << label << " is not set" << endl;
        exit(0);
    }

    label = base_label + "/numOfNodeInElmVelocity";
    if ( !tp.getInspectedValue(label,numOfNodeInElmVelocity)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    label = base_label + "/numOfNodeInElmPressure";
    if ( !tp.getInspectedValue(label,numOfNodeInElmPressure)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    base_label = "/Parameter";
    label = base_label + "/mu";
    if ( !tp.getInspectedValue(label,mu)){
        cout << label << " is not set" << endl;
        exit(0);
    }
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
    label = base_label + "/numOfompthreads";
    if ( !tp.getInspectedValue(label,numOfOmpThreads)){
        cout << label << " is not set" << endl;
        exit(0);
    }
}