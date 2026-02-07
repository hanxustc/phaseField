#include "simPF.h"

void Simulation::readSimParameters(int& argc, char**& argv) {
    int simulation;
    ParameterInt paraInt;
    ParameterMater paraMater;
    InitialCondition init(&paraInt);
    Mat<double> theta = init.theta_init();
    ParameterPF paraPF(theta);
    MPIHandler mpi(argc, argv);
    PhaseField2D pf2d(&paraPF, &paraInt, &paraMater, &init, &mpi);

    ifstream in;
    in.open("simPF.in");
    if (in.fail()) {
        cout << "Cannot open simPF.in!\n";
        return;
    }
    in >> simulation;
    in.close();

    switch (simulation) {
    case 1:
        cout << "Not yet implemented!" << endl;
        break;
    case 2:
        pf2d.fdm();
        break;
    default:
        cout << "Invalid simulation selected! Aborting" << endl;
        break;
    }
    return;
}