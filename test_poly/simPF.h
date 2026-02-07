#ifndef _simPF_h_
#define _simPF_h_

#include "pf.h"

class Simulation {
    PhaseField1D* pf1d;
    PhaseField2D* pf2d;
public:
    ParameterPF* parameterPF;
    ParameterInt* parameterInt;
    ParameterMater* parameterMater;
    InitialCondition* initialCondition;
    MPIHandler* mpiHandler;

    void readSimParameters(int& argc, char**& argv);
};

#endif