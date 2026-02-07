#include "parameters.h"

Parameter::Parameter() {}

ParameterPF::ParameterPF(Mat<double> _theta) : theta(_theta) {}

ParameterInt::ParameterInt() {
    ifstream in;

    in.open("paramIntegration.in");
    if (in.fail()) {
        cout << "Unable to open paramIntegration.in for integration parameters" << endl;
        return;
    }
    in >> dx;
    in >> dy;
    in >> dt;
    in >> nstep;
    in >> nprint;
    in >> Nx;
    in >> Ny;
    in >> err;
    nsmooth = 5000;
    nthetapm = 100;
}

ParameterMater::ParameterMater() {
    ifstream in;

    in.open("paramMaterial.in");
    if (in.fail()) {
        cout << "Unable to open paramMaterial.in for material parameters" << endl;
        return;
    }
    in >> M;
    in >> B0;
    in >> alpha;
    in >> beta;
    in >> u0;
    in >> u1;
    in >> n;
    in >> m;
    in >> epsilon;
    in >> K;
}