#ifndef _pf_h_
#define _pf_h_

#include "initCondition.h"
#include "mpiHandler.h"

class PhaseField {
public:
    PhaseField();
};

class PhaseField1D : public PhaseField {
public:
    PhaseField1D();
};

class PhaseField2D : public PhaseField {
    Mat<double> grad_theta, norm_gradtheta, norm2_gradtheta, Psi;
    Mat<double> theta_smooth, theta_min, theta_max, dtheta, u1_new, alpha_new;
    Mat<double> w, dw, c, dc, B, b, db, dL_theta, dL_gradtheta, df;
    double F1, F2, F3, F = 0;
    int cols_per_process, start_col, end_col;
    void mpiSplit();
    void allocate();
    void smooth();
    void smoothing();
    void gradient();
    void inclination();
    void thetaPlusthetaMinusExtr2G();
    void thetaPlusthetaMinusExtr();
    void newTerm2G();
    void newTerm();
    void weightFunction();
    void gbEnergy();
    void doubleWell();
    void dFunctional();
    void totalEnergy();
public:
    ParameterPF* parameterPF;
    ParameterInt* parameterInt;
    ParameterMater* parameterMater;
    InitialCondition* initialCondition;
    MPIHandler* mpiHandler;
    PhaseField2D(ParameterPF* _parameterPF, ParameterInt* _parameterInt, ParameterMater* _parameterMater,
                 InitialCondition* _initialCondition,
                 MPIHandler* _mpiHandler);
    void fdm();
};

#endif