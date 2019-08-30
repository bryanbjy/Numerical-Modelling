#include <iostream>
#include <string>

//Main class that uses rungekutta to evaluate concentrations x1 and x2 at time t for a range of t values
class CPhosphate {
public:
	CPhosphate();
	CPhosphate(double t_in, double x1_in, double x2_in, double tau_in, double E_in, int step_in, int end_in);
	~CPhosphate();
	double t, x1, x2; //initial conditions of time t, concentration x1 and concentration x2
	double** rkArray; //nested array to contain t array, x1 array and x2 array
	int step, end; //step size of t and upper limit of t array to be evaluated
	double dxdt1(double t, double x1, double x2, double tau_in, double E_in);
	double dxdt2(double t, double x1, double x2, double tau_in, double E_in);
	double noiseG(double x_in);
	double** rungekutta(double t, double x1, double x2, double tau_in, double E_in);
	double** error(CPhosphate rhs);
	void fileWrite(std::string* str1, std::string* str2, std::string* str3);
};


//Inheritance child class for parameter optimisation of tau and reminEff
class optPar : public CPhosphate {
public:
	optPar(double tau_in, double E_in, double tau_opt, double E_opt, int step_in, int end_in);
	double costfunc(double tau_in, double E_in, double tau_opt, double E_opt);
	double* gradDes(double tau_in, double E_in, double tau_opt, double E_opt);
};


//Analysis class for analyzing error simulated from 'real' and 'random' datasets
class Analysis {
public:
	Analysis();
	Analysis(double* t_array, double* x_array, int end_in);
	~Analysis();
	double *timeArray, *xArray; // for storing t and x arrays that will be assigned by the user
	int size; // size of t array
	double* smallestError(double* timeArray, double* xArray, int size_in);
	double* SLR(double* timeArray, double* xArray, int size_in);
	double sqrtErr(double* xArray, double* x_model, int size_in);
	double coDet(double* xArray, double* x_model, int size_in);
};