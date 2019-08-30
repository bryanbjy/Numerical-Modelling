#include "pModel.h"
#include <cmath>
#include <fstream>
#include <string>

using namespace std;


//Stand alone function to read files in a directory to be later analyzed
double** fileRead(int num) {
	int count;
	double **full_array;
	full_array = new double*[3];
	count = 0;
	double t, x1, x2;
	std::string fileName;
	std::string file;

	fstream myFile;
	fileName = new char[5];
	cout << "Which error file do you want to read and analyse?: " << endl;
	cin >> fileName;
	file = "H:\\Visual Studio 2015\\Projects\\pModel\\pModel\\" + fileName + ".txt"; //file is a combination of directory to file folder and the file name the user has specified, must be changed to user's directory

	myFile.open(file, fstream::in);

	for (int i = 0; i < 3; i++) {
		full_array[i] = new double[num];
	}

	if (myFile.fail()) {
		cout << "Error opening file" << endl;
		exit(0);
	}

	std::string line;
	getline(myFile, line);

	//while end of file is not reached, read in elements line by line and put it into arrays
	while (!myFile.eof()) { 
		myFile >> t;
		myFile >> x1;
		myFile >> x2;
		full_array[0][count] = t;
		full_array[1][count] = x1;
		full_array[2][count] = x2;

		count += 1;

	}
	myFile.close();

	for (int i = 0; i < num + 1; i++) {
		cout << full_array[0][i] << "  " << full_array[1][i] << "  " << full_array[2][i] << endl;
	}

	return(full_array);
}

int main() {

	//Building first two 'real' and 'random' datasets by approximating the coupled ODES with rungekutta
	CPhosphate *real, *random;
	double** val, **store;
	cout << "Building first solution" << endl;
	real = new CPhosphate(0, 0, 0, 100, 0.99, 10, 10000);
	cout << "Building second solution" << endl;
	random = new CPhosphate(0, 0, 0, 100, 0.99, 10, 10000);
	cout << "Calculating for error between real and random data" << endl;
	val = real->error(*random);

	//Reads error file for analysis, gives user an option to analyze other files but will not make much sense
	Analysis *x1Analyse, *x2Analyse;
	store = fileRead(1000);
	cout << endl;
	cout << "Displaying data analysis" << endl;
	cout << endl;
	cout << "x1 data results" << endl;
	x1Analyse = new Analysis(store[0], store[1], 1000);
	cout << endl;
	cout << "x2 data results" << endl;
	x2Analyse = new Analysis(store[0], store[2], 1000);

	//Parameter optimisation with different tau and reminEff values, simulates scenario when actual parameters are unknown
	double E_in, tau_in, E_opt, tau_opt;
	tau_in = 99.9;
	E_in = 0.7;
	tau_opt = 100;
	E_opt = 0.99;
	cout << endl;
	cout << "Optimising Parameters of " << "E = " << E_in << ", Tau = " << tau_in << " to get E = " << E_opt << ", Tau = " << tau_opt << endl;
	optPar* optimisation;
	optimisation = new optPar(tau_in, E_in, tau_opt, E_opt, 10, 10000); //recommended alpha = 0.0001 for this setting

}