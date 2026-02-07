#ifndef _parameters_h_
#define _parameters_h_

#include <fstream>
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

constexpr double pi = 3.14159265358979323846;

class Parameter {
public:
    Parameter();
};

class ParameterPF : public Parameter {
public:
    Mat<double> theta;
    ParameterPF(Mat<double> _theta);
};

class ParameterInt : public Parameter {
public:
    double dx, dy;
    double dt;
    long nstep, nprint, nsmooth, nthetapm;
    long Nx, Ny;
    double err;
    ParameterInt();
};

class ParameterMater : public Parameter {
public:
    double M;
    double B0, alpha, beta;
    double u0, u1;
    double epsilon, K;
    int n, m;
    Mat<int> R = { {0, 1},
                   {-1, 0} };
    ParameterMater();
};

#endif