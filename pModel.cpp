#include "pModel.h"
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

using namespace std;

//set global variables to avoid repeated declaration
double V1 = 3 * pow(10, 16); //box volume of shallow ocean
double F_0 = 6 * pow(10, 14); //overturning volume flux
double F_r = 3 * pow(10, 13); //river water flux
double xR = 1.5; //river water phosphate concentration 
double V2 = 1 * pow(10, 18); //box volume of deep ocean
double tau; //residence time of phosphate in the shallow ocean
double reminEff; //remineralisation coefficient for phosphate
double prodtvy; //productivity rate

//template class for turning doubles to strings for file writing
template <typename T>
std::string numString(T num)
{
	std::ostringstream conv;
	conv << num;
	return conv.str();
}

//default constructor
CPhosphate::CPhosphate() {
	t = 0;
	x1 = 0;
	x2 = 0;
	rkArray = 0;
	step = 0;
	end = 0;
}

//first coupled ODE expression
double CPhosphate::dxdt1(double t, double x1, double x2, double tau_in, double E_in) {
	tau = tau_in;
	prodtvy = x1*V1 / tau;
	return (F_r*xR - F_0*x1 + F_0*x2 - prodtvy) / V1;
}

//second coupled ODE expression
double CPhosphate::dxdt2(double t, double x1, double x2, double tau_in, double E_in) {
	tau = tau_in;
	reminEff = E_in;
	prodtvy = x1*V1 / tau;
	return (F_0*x1 - F_0*x2 + reminEff*prodtvy) / V2;
}

//Rungekutta method to approximate for quantities x1 and x2
double** CPhosphate::rungekutta(double t, double x1, double x2, double tau_in, double E_in) {
	double *t_array, *x1_array, *x2_array, **fullarray, k1, h1, k2, h2, k3, h3, k4, h4;
	int count;

	count = 0;
	t_array = new double[end / step];
	x1_array = new double[end / step];
	x2_array = new double[end / step];
	fullarray = new double*[3];

	//while loop for rungekutta approximation until value of t reaches the specified end value
	while (t < end) {
		h1 = step*(dxdt1(t, x1, x2, tau_in, E_in));
		k1 = step*(dxdt2(t, x1, x2, tau_in, E_in));
		h2 = step*(dxdt1(t + 0.5*step, x1 + 0.5*k1, x2 + 0.5*h1, tau_in, E_in));
		k2 = step*(dxdt2(t + 0.5*step, x1 + 0.5*k1, x2 + 0.5*h1, tau_in, E_in));
		h3 = step*(dxdt1(t + 0.5*step, x1 + 0.5*k2, x2 + 0.5*h2, tau_in, E_in));
		k3 = step*(dxdt2(t + 0.5*step, x1 + 0.5*k2, x2 + 0.5*h2, tau_in, E_in));
		h4 = step*(dxdt1(t + step, x1 + k3, x2 + h3, tau_in, E_in));
		k4 = step*(dxdt2(t + step, x1 + k3, x2 + h3, tau_in, E_in));

		x1 = x1 + (1. / 6.)*(h1 + 2.*h2 + 2.*h3 + h4);
		x2 = x2 + (1. / 6.)*(k1 + 2.*k2 + 2.*k3 + k4);
		t = t + step;
		x1_array[count] = x1;
		x2_array[count] = x2;
		t_array[count] = t;
		count += 1;

	}

	//assigning t, x1 and x2 array into a nested array called fullarray so that it can be returned in the member function
	fullarray[0] = t_array;
	fullarray[1] = x1_array;
	fullarray[2] = x2_array;

	return(fullarray);
}

//method to generate "random" noise in a dataset for comparison to "real" data
double CPhosphate::noiseG(double x_in) {
	double f;
	int num;

	f = (double)rand() / RAND_MAX;

	num = rand() % 2;

	//random number multiplied by a simple scaling to simulate noise 
	if (num == 0) {
		x_in += f*(x_in / 10);
	}
	else {
		x_in -= f*(x_in / 10);
	}

	return x_in;
}

//method to write files
void CPhosphate::fileWrite(std::string* str1, std::string* str2, std::string* str3) {
	std::string fileName; 
	std::string file;

	fstream myFile;
	fileName = new char[5];
	cout << "Writing into file; what do you want your file name to be:" << endl; //user input for name of file
	cin >> fileName;
	file = "H:\\Visual Studio 2015\\Projects\\pModel\\pModel\\" + fileName + ".txt"; //file is a combination of directory to file folder and the file name the user has specified, must be changed to user's directory


	myFile.open(file, fstream::out);


	if (myFile.fail()) {
		cout << "Error opening file" << endl;
		exit(0);
	}

	//Header of text file
	myFile << "Time(t/yr)" << "  ";
	myFile << "Concentration(x1/mmol)" << "  ";
	myFile << "Concentration(x2/mmol)" << endl;

	//writes each element from string arrays str1, str2 and str3 into a file seperated by double spaces
	for (int i = 0; i < end / step; i++) {
		myFile << str1[i] << "  ";
		myFile << str2[i] << "  ";
		myFile << str3[i] << endl;
		cout << str1[i] << "  " << str2[i] << "  " << str3[i] << endl;
	}

	myFile.close();
}

//Constructor which automatically calls rungekutta method to approximate x1 and x2 for time t with given number of steps and end of iteration
CPhosphate::CPhosphate(double t_in, double x1_in, double x2_in, double tau_in, double E_in, int step_in, int end_in) {
	srand(time(NULL));
	std::string randStr;

	t = t_in; 
	x1 = x1_in;
	x2 = x2_in;
	step = step_in;
	end = end_in;
	cout << "Is this real or random data (real/rand)? " << endl; //ask user if the dataset should be the model or the experimental dataset
	cin >> randStr;

	rkArray = new double*[3];
	std::string* stringlist1; 
	std::string* stringlist2;
	std::string* stringlist3;

	stringlist1 = new std::string[end / step]; 
	stringlist2 = new std::string[end / step];
	stringlist3 = new std::string[end / step];

	cout << "dxdt1: " << dxdt1(t, x1, x2, tau_in, E_in) << endl;
	cout << "dxdt2: " << dxdt2(t, x1, x2, tau_in, E_in) << endl;

	for (int j = 0; j < 3; j++) {
		rkArray[j] = new double[end / step];
	}

	rkArray = rungekutta(t, x1, x2, tau_in, E_in); //assigning nested array from rungekutta evaluation to class attribute rkArray so datapoints could be stored and passed

	for (int j = 0; j < 3; j++) {
		for (int i = 0; i <end / step; i++) {
			cout << "Runge Kutta " << j << "; step " << i + 1 << ":" << rkArray[j][i] << endl;
		}
	}

	if (randStr == "rand") {
		for (int j = 1; j < 3; j++) { //don't perform noise on time array
			for (int i = 0; i < end / step; i++) {
				rkArray[j][i] = noiseG(rkArray[j][i]);
			}
		}
	}


	for (int i = 0; i < end / step; i++) {

		stringlist1[i] = numString(rkArray[0][i]);
		stringlist2[i] = numString(rkArray[1][i]);
		stringlist3[i] = numString(rkArray[2][i]);

	}

	fileWrite(stringlist1, stringlist2, stringlist3); //fileWrite function to write files of data points stored within rkArray

}

//destructor to destroy main storing array with t, x1, and x2 arrays nested inside
CPhosphate::~CPhosphate() {
	for (int j = 0; j < 3; j++) {
		delete[] rkArray[j];
	}
	delete[] rkArray;
}

//Compares to objects from class CPhosphate and gives the resulting difference in a nested array
double** CPhosphate::error(CPhosphate rhs) {
	double** temp;
	temp = new double*[3];
	std::string* stringlist1;
	std::string* stringlist2;
	std::string* stringlist3;

	stringlist1 = new std::string[end / step];
	stringlist2 = new std::string[end / step];
	stringlist3 = new std::string[end / step];


	for (int j = 0; j < 3; j++) {
		temp[j] = new double[end / step];
	}

	for (int i = 0; i < end / step; i++) {
		temp[0][i] = rhs.rkArray[0][i]; //setting time array of CPhosphate object as temp object's time array as error evaluation does not concern t
	}

	for (int j = 1; j < 3; j++) {
		for (int i = 0; i < end / step; i++) {
			temp[j][i] = abs(rkArray[j][i] - rhs.rkArray[j][i]); //resulting difference of each data point in x1 and x2 array from model dataset and experimental dataset
		}
	}


	for (int i = 0; i < end / step; i++) {

		stringlist1[i] = numString(temp[0][i]);
		stringlist2[i] = numString(temp[1][i]);
		stringlist3[i] = numString(temp[2][i]);

	}

	fileWrite(stringlist1, stringlist2, stringlist3);

	return(temp);
}


//Cost function to minimise
double optPar::costfunc(double tau_in, double E_in, double tau_opt, double E_opt) {
	double* xObs1, *xReal1, *xObs2, *xReal2;
	double cost;

	cost = 0;

	xObs1 = new double[end / step]; //x1 value from initial guess parameter
	xReal1 = new double[end / step]; //x1 value from actual parameter we want to eventually optimize
	xObs2 = new double[end / step]; //x2 value from initial guess parameter
	xReal2 = new double[end / step]; //x2 value from actual parameter we want to eventually optimize

	xObs1 = rungekutta(0, 0, 0, tau_in, E_in)[1];
	xObs2 = rungekutta(0, 0, 0, tau_in, E_in)[2];
	xReal1 = rungekutta(0, 0, 0, tau_opt, E_opt)[1];
	xReal2 = rungekutta(0, 0, 0, tau_opt, E_opt)[2];

	for (int i = 0; i < end / step; i++) {
		cost += pow((xObs1[i] - xReal1[i]), 2) + pow((xObs2[i] - xReal2[i]), 2); //least square expression as cost function
	}
	cout << cost << endl;
	return(cost);
}

//Optimisation method: Gradient Descent - has been tested with the Rosenbrock function -> See Appendix
double* optPar::gradDes(double tau_in, double E_in, double tau_opt, double E_opt) {
	int counts;
	double val_t, val_E, rel_imp, dt, dE, dJdt, dJdE, alpha, J_new, J_old;
	double* parContainer;
	val_t = tau_in; //initial guess tau parameter
	val_E = E_in; //initial guess remineralisation coefficient parameter
	tau = tau_opt; //actual tau parameter we want to eventually optimize
	reminEff = E_opt; //actual remineralisation coefficient parameter we want to eventually optimize
	rel_imp = 1; 
	dt = tau_in / 10000; //scale four magnitudes down
	dE = E_in / 10000; //scale four magnitudes down
	counts = 0; //to count number of iteration before reaching tolerance

	cout << "Please key in the tuning parameter alpha: " << endl; //needs manual tuning
	cin >> alpha;

	//while loop to evaluate cost function with gradient descent method until relative improvement is smaller than 10^-6
	while (rel_imp > pow(10, -6)) { 
		counts += 1;
		J_old = costfunc(val_t, val_E, tau, reminEff); 
		dJdt = (costfunc(val_t + dt, val_E, tau, reminEff) - costfunc(val_t - dt, val_E, tau, reminEff)) / (2 * dt); //approximating gradient of cost function with respect to tau parameter with central differencing
		dJdE = (costfunc(val_t, val_E + dE, tau, reminEff) - costfunc(val_t, val_E - dE, tau, reminEff)) / (2 * dE); //approximating gradient of cost function with respect to E parameter with central differencing
		val_t = val_t - alpha*dJdt; //new tau parameter to evaluate on next cost function call and next iteration in the while loop 
		val_E = val_E - alpha*dJdE; //new E parameter to evaluate on next cost function call and next iteration in the while loop
		J_new = costfunc(val_t, val_E, tau, reminEff); //resultant cost function from new tau and E parameters
		rel_imp = abs(J_old - J_new) / ((J_old + J_new) / 2); //calculating relative improvement that is used as the condition of this while loop
		cout << "count " << counts << ": " << rel_imp << endl;
	}

	parContainer = new double[3];
	parContainer[0] = val_t;
	parContainer[1] = val_E;
	parContainer[2] = J_new;

	cout << "Optimised tau parameter: " << val_t << endl;
	cout << "Optimised remineralisation coefficient parameter: " << val_E << endl;
	cout << "Terminated cost function value: " << J_new << endl;

	return(parContainer);
}

//constructor calls for gradient descent function
optPar::optPar(double tau_in, double E_in, double tau_opt, double E_opt, int step_in, int end_in) {
	double* J;
	t = 0;
	x1 = 0;
	x2 = 0;
	step = step_in;
	end = end_in;
	J = new double[3];
	J = gradDes(tau_in, E_in, tau_opt, E_opt);
}


Analysis::Analysis() {
	timeArray = 0;
	xArray = 0;
	size = 0;
}

Analysis::~Analysis() {
	delete[] timeArray;
	delete[] xArray;
}

//method to find minimum error with O(n) time complexity
double* Analysis::smallestError(double* timeArray, double* xArray, int size_in) {
	double first, second;
	double* errorArray;

	first = xArray[0]; //first element in array that contains the smallest value 
	second = xArray[0]; //second element in array that contains the second smallest value 
	errorArray = new double[2];

	for (int i = 0; i < size_in; i++) {
		if (xArray[i] < first) { //when x element from xArray is smaller than first, reassign first value as second value and the x element as the first value
			second = first;
			first = xArray[i];
			errorArray[0] = timeArray[i];
			errorArray[1] = xArray[i];
		}

		else if (xArray[i] < second && xArray[i] > first) { //when x element from xArray is smaller than second but bigger than first, reassign second value as x element
			second = xArray[i];
		}
	}

	cout << "SmallestErr - Minimum error in time: " << errorArray[0] << endl;
	cout << "SmallestErr - Minimum error of x-value: " << errorArray[1] << endl;

	return(errorArray);

}

//Simple Linear Regression
double* Analysis::SLR(double* timeArray, double* xArray, int size_in) {
	double meanT, meanTsqrt, meanTX, meanX, m, b;
	double* x_model;

	meanT = 0; //the mean of timeArray
	meanTsqrt = 0; //the mean of timeArray squared
	meanTX = 0; //the mean of time Array multiplied by xArray
	meanX = 0; //the mean of xArray
	x_model = new double[size_in];

	for (int i = 0; i < size_in; i++) {
		meanT += timeArray[i]; 
		meanX += xArray[i];
		meanTsqrt += timeArray[i] * timeArray[i];
		meanTX += timeArray[i] * xArray[i];
	}

	meanT = meanT / size_in;
	meanX = meanX / size_in;
	meanTsqrt = meanTsqrt / size_in;
	meanTX = meanTX / size_in;

	m = (meanT * meanX - meanTX) / (meanT*meanT - meanTsqrt);
	b = meanX - m * meanT;

	cout << "Simple Linear Regression model" << endl;
	cout << "m value: " << m << endl;
	cout << "b value: " << b << endl;
	for (int i = 0; i < size_in; i++) { //calculating for resultant x values from model by Simple Linear Regression
		x_model[i] = m*timeArray[i] + b;
	}

	return(x_model);
}

//R squared error calculation for line of best fit
double Analysis::sqrtErr(double* xArray, double* x_model, int size_in) {
	double R;
	R = 0;

	for (int i = 0; i < size_in; i++) {
		R += (x_model[i] - xArray[i])*(x_model[i] - xArray[i]);
	}

	return R;
}

//Coefficient of determination to quantify the best fit
double Analysis::coDet(double* xArray, double* x_model, int size_in) {
	double* x_meanl;
	double sqrt_regr, sqrt_xmeanl;
	double meanX;
	double D;

	meanX = 0;
	x_meanl = new double[size_in];

	for (int i = 0; i < size_in; i++) {
		meanX += xArray[i];
	}

	meanX = meanX / size_in;

	for (int i = 0; i < size_in; i++) {
		x_meanl[i] = meanX;
	}

	sqrt_regr = sqrtErr(xArray, x_model, size_in); //R squared error between x SLR model and x dataset 
	sqrt_xmeanl = sqrtErr(xArray, x_meanl, size_in); //R squared error between mean of x dataset and x dataset 
	D = 1 - (sqrt_regr / sqrt_xmeanl);

	cout << "R squared error: " << sqrt_regr << endl;
	cout << "Coefficient of Determination: " << D << endl;

	return D;
}

//Constructor calls for all analysis member functions
Analysis::Analysis(double* t_array, double* x_array, int size_in) {
	double* minErr, *slrModel;
	double Rsqrt, Det;
	timeArray = t_array;
	xArray = x_array;
	size = size_in;

	minErr = new double[2];
	minErr = smallestError(timeArray, xArray, size); 
	slrModel = new double[size]; //array slrModel is initialised to contain the resulting x values from the SLR simple linear regression function
	slrModel = SLR(timeArray, xArray, size); 
	Rsqrt = sqrtErr(xArray, slrModel, size);
	Det = coDet(xArray, slrModel, size);

}
