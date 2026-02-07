#include "initCondition.h"

InitialCondition::InitialCondition(ParameterInt* _parameterInt) : parameterInt(_parameterInt) {
    read_theta();
    vec uniqueElements = unique(vectorise(theta));
    theta_num = uniqueElements.t();
    ngrain = theta_num.size();
}

void InitialCondition::read_theta() {

    ifstream in;
    in.open("initialCondition.in");
    if (in.fail()) {
        cout << "Unable to open initialCondition.in for initial conditions" << endl;
        return;
    }
    
    ifstream in_theta;
    string ic;
    in >> ic;
    if (ic == "poly") {
        in_theta.open("poly_matrix.in");
        if (in_theta.fail()) {
            cout << "Unable to read theta matrix" << endl;
            return;
        }
        mat poly_matrix;
        theta.load("poly_matrix.in");
    }
    else if (ic == "triple") {
        in_theta.open("triple_matrix.in");
        if (in_theta.fail()) {
            cout << "Unable to read theta matrix" << endl;
            return;
        }
        mat triple_matrix;
        theta.load("triple_matrix.in");
    }
    else if (ic == "Wulff" || ic == "growth") {
        in_theta.open("circle_matrix.in");
        if (in_theta.fail()) {
            cout << "Unable to read theta matrix" << endl;
            return;
        }
        mat circle_matrix;
        theta.load("circle_matrix.in");
    }
    else if (ic == "energy") {
        in_theta.open("energy_matrix.in");
        if (in_theta.fail()) {
            cout << "Unable to read theta matrix" << endl;
            return;
        }
        mat energy_matrix;
        theta.load("energy_matrix.in");
    }
    else {
        cout << "Invalid initial condition" << endl;
        return;
    }
}

Mat<double> InitialCondition::theta_init() {
    read_theta();
    return theta;
}