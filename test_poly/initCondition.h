#ifndef _initCondition_h_
#define _initCondition_h_

#include "parameters.h"

class InitialCondition {
public:
	ParameterInt* parameterInt;
	InitialCondition(ParameterInt* _parameterInt);
	void read_theta();

	int ngrain;
	Mat<double> theta_init();
	Mat<double> theta;
	Mat<double> theta_num;
};

#endif