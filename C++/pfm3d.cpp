#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <new>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <cstdarg>
#include <string>
using namespace std;

double *NewRho;
double **eta, **circleProp, **BFE_partialEta, **NewEta, **SwitchEta;
double *rho, *BFE_partialRho, *Dif, *h, *SwitchRho, *rhoDif;
double RaND, etaSqSum, etaCubeSum, rhoSq, rhoCube, rhoQuad, etaSqZero, etaSqOne;
double phi, etaCubeZero, etaCubeOne, sumrho, sumeta, sumrhoInit = 0, sumetaInit = 0, intvrho;
double t, finalt, S, SpR, Sdif, evat, evatf, tempr, intvs, intvf;
int Size, N_size, countt, counttf;
int pixelRadius, rad, sizebv;
int **N;

int i0, i1, i2, i3, i4, i5, i6, i7;

double T, seed = 0;

double A = 12, B = 1, GD = 4.1, GE = 3.75, GM = 10, SD = 45, SE = 10, VD = 0.08, prev = 0;
double MA = 474, MB = 84, MGD = 24, MGE = 165, MGM = 355, MSD = 442, MSE = 44, MVD = 13;

double uba = MA, ubb = MB, ubgd = MGD, ubge = MGE, ubgm = MGM, ubsd = MSD, ubse = MSE, ubvd = MVD; // upper bounds for constants
double lba = 0, lbb = 0, lbgd = 0, lbge = 0, lbgm = 0, lbsd = 0, lbse = 0, lbvd = 0;			   // lower bounds for constants

//Arrays of parameter values
double params[] = {A, B, GD, GE, GM, SD, SE, VD};					  // Array of values used in the calculation
double MAX[] = {MA, MB, MGD, MGE, MGM, MSD, MSE, MVD};				  // Array of the maximum for each parameter
double upperbound[] = {uba, ubb, ubgd, ubge, ubgm, ubsd, ubse, ubvd}; // Array to track parameter upperbound

double lowerbound[] = {lba, lbb, lbgd, lbge, lbgm, lbsd, lbse, lbvd};   // Array to track parameter lowerbound
string mapback_dict[] = {"A", "B", "GD", "GE", "GM", "SD", "SE", "VD"}; // Array used check the current parameter for reassigning values

int current; // the current running combination
//grouping relevant data together
struct variables
{
	int xlength, ylength, zlength;
	int noOfparticles;
	double rhoVap, etaVap, rhoSolid, etaSolid;
} dim;

struct CircleInfo
{
	int midX, midY, midZ, particleNo;
	double value, radius;
} circle;

typedef struct constants
{
	double surfdif, gbdif, voldif, vapdif;
	double surfenergy, gbenergy, gbmobility;
	double randR, randE, deltaT;
	double A, B, Cgbe;
} constants;
constants con;

//memory allocation for 1D
double *store1D(int maxArray)
{
	double *k;
	k = new double[maxArray];

	if (k == NULL)
	{
		cout << "Failed to allocate 2D memory";
		exit(EXIT_FAILURE);
	}
	return k;
}

//memory allocation for 2D
double **store2D(int maxArray1, int maxArray2)
{
	double **m;
	m = new double *[maxArray1];

	for (int a = 0; a < maxArray1; a++)
	{
		m[a] = new double[maxArray2];
	}

	if (m == NULL)
	{
		cout << "Failed to allocate 2D memory";
		exit(EXIT_FAILURE);
	}
	return m;
}

int **store2Dint(int maxArray1, int maxArray2)
{
	int **m;
	m = new int *[maxArray1];

	for (int a = 0; a < maxArray1; a++)
	{
		m[a] = new int[maxArray2];
	}

	if (m == NULL)
	{
		cout << "Failed to allocate 2D memory";
		exit(EXIT_FAILURE);
	}
	return m;
}

//clear 1D array
void clear1D(double *j1)
{
	delete[] j1;
	return;
}

//clear 2D array
void clear2D(double **j2, int a)
{
	for (int i = 0; i < a; i++)
	{
		delete[] j2[i];
	}
	delete[] j2;
	return;
}

void clear2Dint(int **j2, int a)
{
	for (int i = 0; i < a; i++)
	{
		delete[] j2[i];
	}
	delete[] j2;
	return;
}

/*Draw circle to assign the field variables in the circle (rho and eta for solid) and the variables outside the circle in the vapor stage */
void CreateParticles(double *g, struct variables dim, struct CircleInfo circle)
{
	int x, y, z, i;

	for (z = 0; z <= dim.zlength; z++)

	{
		for (y = 0; y <= dim.ylength; y++)
		{
			for (x = 0; x < dim.xlength; x++)
			{
				if ((x - circle.midX) * (x - circle.midX) + (y - circle.midY) * (y - circle.midY) + (z - circle.midZ) * (z - circle.midZ) <= circle.radius * circle.radius)
				{
					i = (dim.xlength * dim.ylength * z) + (dim.xlength * y) + x;

					g[i] = circle.value;
				}
			}
		}
	}

	return;
}

/*print results to a file */

void fileprintfull(double t /*like the actual time t*/, double timet /*time you want*/, double deltatime, double *R, double **E, string namestart, string nameend /*include .txt*/, double evalt, int Size)

{
	ofstream myfile;

	double f1 = timet - deltatime;
	double f2 = timet + deltatime;

	if ((t > f1) && (t < f2))
	{

		ostringstream oss;
		string var;
		char *name;
		oss << namestart << evalt << nameend;

		var = oss.str();

		name = new char[var.length()];

		strcpy(name, var.c_str());
		ofstream myfile;

		myfile.open(name);

		for (int z = 0; z < Size; z++)
		{
			myfile << R[z] << "\t" << E[0][z] << "\t" << E[1][z] << "\n";
		}

		myfile.close();
	}
}

void fileprintsum /*matlab form*/ (double t /*like the actual time t*/, double timet /*time you want*/, double deltatime, double *R, double **E, string namestart, string nameend /*include .txt*/, double evalt, int Size)

{

	double f1 = timet - deltatime;
	double f2 = timet + deltatime;

	if ((t > f1) && (t < f2))
	{
		double sum;
		double sump;

		ostringstream oss;
		string var;
		char *name;
		oss << namestart << evalt << nameend;

		var = oss.str();

		name = new char[var.length()];

		strcpy(name, var.c_str());
		ofstream myfile;

		myfile.open(name);

		for (int z = 0; z < Size; z++)
		{

			sum = R[z] * (E[0][z] + E[1][z]);

			myfile << sum << "\n";
		}

		myfile.close();
	}
}

void fileprintrho /*matlab form*/ (double t /*like the actual time t*/, double timet /*time you want*/, double deltatime, double *R, string namestart, string nameend /*include .txt*/, double evalt, int Size)

{

	double f1 = timet - deltatime;
	double f2 = timet + deltatime;

	if ((t > f1) && (t < f2))
	{
		double sum;
		double sump;

		ostringstream oss;
		string var;
		char *name;
		oss << namestart << evalt << nameend;

		var = oss.str();

		name = new char[var.length()];

		strcpy(name, var.c_str());
		ofstream myfile;

		myfile.open(name);

		for (int z = 0; z < Size; z++)
		{
			myfile << R[z] << "\n";
		}

		myfile.close();
	}
}
// this function has no importance in this program
// but check first before you purge
void printUsage()
{
	char usage[] = {"To use the program:\nSpecify the number of constants followed by their values and simulation time\nIn the\
 following order:\n#number, A, B, gbdif, gbenergy, gbmobility, surfdif, surfenergy, voldif, finalt \0"};
	cout << usage << endl;
	exit(-1);
}

double inRange(double low, double high, double seed)
{
	// this function is magic that makes picking random numbers easy
	// low: smallest number in the range
	// high: greatest number in range
	// seed: to seed the srand function below => Do I even need to tell that? hope not! :0
	srand(seed);
	return (low + (high - low) * (double)rand() / RAND_MAX);
}

void resetParams()
{
	// Does what its name says
	// so deal with it

	con.A = 12.0;
	con.B = 1.0;
	con.gbdif = 4.1;
	con.gbenergy = 3.75;
	con.gbmobility = 10.0;
	con.surfdif = 45.0;
	con.surfenergy = 10.0;
	con.voldif = 0.08;
}
void resetRho(double *rho, double *newrho, double **eta, int Size)
{
	// rho: pointer to the density array
	// newrho: pointer to the new calculated density array
	// Size: size of the rho array
	// sets the values of rho and eta to default values and newrho to zero

	srand(0);
	static int count = 0;
	sumrhoInit = sumetaInit = 0;
	//start by setting the entire sample area to the vapor phase
	memset(newrho, 0, Size);
	for (int i = 0; i < Size; i++)
	{
		double RaND = ((double)rand() / RAND_MAX);
		rho[i] = dim.rhoVap + (0.5 - RaND) * con.randR;
		//density field variable
		for (int g = 0; g < dim.noOfparticles; g++)
		{
			eta[g][i] = dim.etaVap + (0.5 - RaND) * con.randE; //orientation parameter
		}
	}

	rad = pixelRadius;

	circleProp[0][0] = (double)rad;
	circleProp[1][0] = (double)rad / 2;
	circleProp[0][1] = 15;
	circleProp[1][1] = 31;
	circleProp[0][2] = dim.ylength / 2;
	circleProp[1][2] = dim.ylength / 2;
	circleProp[0][3] = dim.zlength / 2;
	circleProp[1][3] = dim.zlength / 2;

	//filling the circle with the right field variable value for rho and eta
	for (int g = 0; g < dim.noOfparticles; g++)
	{
		circle.radius = circleProp[g][0];
		circle.midX = circleProp[g][1];
		circle.midY = circleProp[g][2];
		circle.midZ = circleProp[g][3];
		circle.value = dim.rhoSolid;
		circle.particleNo = g;
		CreateParticles(rho, dim, circle);
		circle.value = dim.etaSolid;
		CreateParticles(eta[circle.particleNo], dim, circle);
	}
}
// these global variables are very important
// delete at your own headache
double RHO = 4700;
double TOL = 0.001;
void reassign(constants *cons, int count, ...)
{
	// This method job is to copy value into the variables used in the calculation
	// con: pointer to the struct of the constants
	// count: number of variable arguments passed
	// ... : variable arguments
	va_list args;
	va_start(args, count);

	for (int i = 0; i < count; i++)
	{
		int index = va_arg(args, int);

		if (index == 0)
			cons->A = params[index];

		else if (index == 1)
			cons->B = params[index];

		else if (index == 2)
			cons->gbdif = params[index];

		else if (index == 3)
			cons->gbenergy = params[index];

		else if (index == 4)
			cons->gbmobility = params[index];

		else if (index == 5)
			cons->surfdif = params[index];

		else if (index == 6)
			cons->surfenergy = params[index];

		else if (index == 7)
			cons->voldif = params[index];
	}
	va_end(args);
}
bool modify(const double *sumrho, double *param, double *lowerbound, double *upperbound, const double *maximum, double *t, const int master = 0)
{
	// sumrho: pointer to the density value
	// param: pointer to a parameter
	// lowerbound: pointer to the current lower bound of the parameter
	// upperbound: pointer to the upper bound of the parameter
	// maximum: pointer to the maximum value of the parameter
	// t: pointer to the current time
	// master: determines which method call resets the time because resetting the time is critical for the simulation

	// returns a boolean to either restart or continue with current values of your constants
	static int count = 0; /* this variable is used to prevent an infinite cycle where the lower bound and upper bound are very close
	 and all the values in the range produces an infeasible density value*/
	static bool isWritable = true; // this variable is used to prevent multiple writes to file after the maximum has been found.

	if (*upperbound <= 0 || *upperbound > *maximum)/* checks if the upperbound is less than or equal to zero or greater than the
	maximum. if so sets the upperbound to maximum. Don't remember why this was necessary!*/
		*upperbound = *maximum;										

	if (isnan((*sumrho)) || *sumrho < 0 || *sumrho > RHO)
	{/* checks if the density is infeasible */

		if (*param < *upperbound)
		{// sets the upperbound to the parameter value that produced an infeasible density
			*upperbound = *param;
		}

		cout << "time : " << *t << endl;
		cout << "expecting nan param:" << *param << " rho: " << *sumrho << endl;
		cout << "lowerbound: " << *lowerbound << " upperbound: " << *upperbound << endl;

		*param = inRange(*lowerbound, *upperbound, *t); //randomly picks a value between the upper and lower bounds

		if (count > 5)
		{// checks if count is greater than 5
			*upperbound -= 0.1; // decrements the upper bound by 0.1
			*lowerbound = 0; // sets the lower bound to 0
		}

		if (master)
		{// checks if the current call is the master call
			*t = 0; // sets time to 0 to restart simulation
			if (count > 5) // checks if count is greater than 5
				count = 0; // sets count to 0
			else
				count++; // increments count by one
		}

		return true; // returns true to signal that the simulation needs to be restarted
	}
	else
	{

		if (((*param - *lowerbound) < TOL && master) && isWritable)
		{/* checks if the parameter value is close to its lower bound and that the current call is the master call
			and it is okay to write to file*/
			*t = 0; // sets time to 0

			// writes parameter values to file
			ofstream fout("Combinations.txt", ios::out | ios::app);
			fout << "con.A = " << con.A << ";\n"
				 << "con.B = " << con.B << ";\n"
				 << "con.gbdif = " << con.gbdif << ";\n"
				 << "con.gbenergy = " << con.gbenergy << endl;
			fout << "con.gbmobility = " << con.gbmobility << ";\n"
				 << "con.surfdif = " << con.surfdif << ";\n"
				 << "con.surfenergy = " << con.surfenergy << endl;

			fout << "con.voldif = " << con.voldif << endl;
			fout << "rho = " << *sumrho << endl;
			fout << "=========================" << endl;

			fout.flush();
			fout.close();
			count = current = 0;

			isWritable = false;
			return false;
		}

		if (*param > *lowerbound && *t > T)
		{// checks if the current parameter value is greater than the lower bound and that the time is greater than seed time
			*lowerbound = *param; // sets lowerbound to the current parameter value
			*param = inRange(*lowerbound, *upperbound, *t); // sets the parameter value to a random number between upper and lower bounds
			if (*param < 0 || *param > *maximum)/* checks if the parameter value is less than zero or greater than the maximum*/
				*param = 0; // sets parameter to zero

			if (master)
			{// checks if the current call is the master call 
				*t = 0; // sets time to 0
				if (count < 0) //checks if count is less than zero
					count = 0; // sets count to 0
				else
					count--; // decrement count by 1
			}
			cout << "time : " << *t << endl;
			cout << "expecting non nan param:" << *param << " rho: " << *sumrho << endl;
			cout << "lowerbound: " << *lowerbound << " upperbound: " << *upperbound << endl;

			return true;
		}

		return false;
	}
}
template <typename T>
void clear_1D(T *ptr)
{
	delete[] ptr;
}
template <typename T>
void clear_2D(T **ptr, int size)
{

	for (int i = 0; i < size; i++)
		free(ptr[i]);

	free(ptr);
}
int main(int argc, char *const argv[])
{
	ofstream myfile, myfile2, myfile1;
	myfile.open("fullT0.txt");
	myfile1.open("3p0.txt");
	myfile2.open("sum rho and eta.txt");

	bool writeonce = true; // making sure that combination heading is written once

	//assigning values

	t = 0.0;
	finalt = 1.5001;
	N_size = 6;

	dim.etaSolid = 1.0;
	dim.etaVap = 0.0;
	dim.noOfparticles = 2;
	dim.rhoSolid = 0.9998;
	dim.rhoVap = 0.000000089;
	dim.xlength = 40;
	dim.ylength = 30;
	dim.zlength = 30;

	con.randE = 0.000000;
	con.randR = 0.0000001;
	con.deltaT = 0.0001;
	con.vapdif = 0.0;
	con.Cgbe = 7.0;
// these logic sections extracts the commandline arguments
	if (argc == 5)
	{
		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
	}
	else if (argc == 6)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
	}
	else if (argc == 7)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
		i2 = strtod(argv[6], NULL);
	}
	else if (argc == 8)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
		i2 = strtod(argv[6], NULL);
		i3 = strtod(argv[7], NULL);
	}
	else if (argc == 9)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
		i2 = strtod(argv[6], NULL);
		i3 = strtod(argv[7], NULL);
		i4 = strtod(argv[8], NULL);
	}
	else if (argc == 10)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
		i2 = strtod(argv[6], NULL);
		i3 = strtod(argv[7], NULL);
		i4 = strtod(argv[8], NULL);
		i5 = strtod(argv[9], NULL);
	}
	else if (argc == 11)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
		i2 = strtod(argv[6], NULL);
		i3 = strtod(argv[7], NULL);
		i4 = strtod(argv[8], NULL);
		i5 = strtod(argv[9], NULL);
		i6 = strtod(argv[10], NULL);
	}
	else if (argc == 12)
	{

		finalt = strtod(argv[1], NULL);
		seed = strtod(argv[2], NULL);
		current = strtod(argv[3], NULL);
		i0 = strtod(argv[4], NULL);
		i1 = strtod(argv[5], NULL);
		i2 = strtod(argv[6], NULL);
		i3 = strtod(argv[7], NULL);
		i4 = strtod(argv[8], NULL);
		i5 = strtod(argv[9], NULL);
		i6 = strtod(argv[10], NULL);
		i7 = strtod(argv[11], NULL);
	}
	else
	{ // print usage message
		printf("Usage:\n./program_name [simularion_time  bounding_time sampling_size sample_indices..]\nSample_indices =[0,7]\n");
		printf("sample index maps to this array: [A, B, gbdif, gbenergy, gbmobility, surfdif, surfenergy, voldif]\n");
		exit(-1);
	}
	T = seed; // The time to starting constraining the lowerbound. This ensures that the lowerbound is not constrained before a nan value is found
	// for simulations that fail higher time, this value should be at least half the simulation time to get valid constants that won't
	// fail at higher simulation time. Contants found ara usually good for the same amount of simulation time they were found with.

	con.A = A;			 // 12.0;
	con.B = B;			 //1.0;
	con.gbdif = GD;		 //4.1;
	con.gbenergy = GE;   //3.75;
	con.gbmobility = GM; //10.0;
	con.surfdif = SD;	//45.0;
	con.surfenergy = SE; // 10.0;
	con.voldif = VD;	 // 0.08;

	pixelRadius = 10;
	Size = dim.xlength * dim.ylength * dim.zlength;

	//allocating memory
	rho = store1D(Size * 10);
	BFE_partialRho = store1D(Size);
	Dif = store1D(Size);
	NewRho = store1D(Size);
	h = store1D(Size);
	SwitchRho = store1D(Size);
	rhoDif = store1D(Size);

	eta = store2D(dim.noOfparticles, Size);
	circleProp = store2D(dim.noOfparticles, N_size);
	BFE_partialEta = store2D(dim.noOfparticles, Size);
	NewEta = store2D(dim.noOfparticles, Size);
	SwitchEta = store2D(dim.noOfparticles, Size);

	N = store2Dint(N_size, Size);

R:// Label for goto statement

	//start by setting the entire sample area to the vapor phase
	resetRho(rho, NewRho, eta, Size);

	//cout<<"check\n";
	//check initial circle microstructure

	for (int z = 0; z < dim.zlength; z++)
	{
		for (int y = 0; y < dim.ylength; y++)
		{
			for (int x = 0; x < dim.xlength; x++)
			{

				int l = (dim.xlength * dim.ylength * z) + (y * dim.xlength) + x;

				myfile << rho[l] << "\t" << eta[0][l] << "\t" << eta[1][l] << "\n";
			}
		}
	}

	for (int z = 0; z < dim.zlength; z++)
	{
		for (int y = 0; y < dim.ylength; y++)
		{
			for (int x = 0; x < dim.xlength; x++)
			{

				int l = (dim.xlength * dim.ylength * z) + (y * dim.xlength) + x;

				myfile1 << rho[l] << "\n";
			}
		}
	}

	for (int i = 0; i < Size; i++)
	{
		sumrhoInit = sumrhoInit + rho[i];

		for (int g = 0; g < dim.noOfparticles; g++)
		{
			sumetaInit = sumetaInit + eta[g][i];
		}
	}

	//Creating nearest neigbours for numerical differenciation using partial derivatives
	for (int z = 0; z <= dim.zlength; z++)
	{
		for (int y = 0; y <= dim.ylength; y++)
		{
			for (int x = 0; x < dim.xlength; x++)
			{
				int k;

				k = (dim.xlength * dim.ylength * z) + (y * dim.xlength) + x;

				N[0][k] = ((x + dim.xlength + 1) % dim.xlength) + ((y + dim.ylength) % dim.ylength) * dim.xlength + (dim.ylength * dim.xlength * z);								 //f(x+dx,y,z) dx=dy=dz=di=1
				N[1][k] = ((x + dim.xlength - 1) % dim.xlength) + ((y + dim.ylength) % dim.ylength) * dim.xlength + (dim.ylength * dim.xlength * z);								 //f(x-dx,y,z)
				N[2][k] = ((x + dim.xlength) % dim.xlength) + ((y + dim.ylength + 1) % dim.ylength) * dim.xlength + (dim.ylength * dim.xlength * z);								 //f(x,y+dy,z)
				N[3][k] = ((x + dim.xlength) % dim.xlength) + ((y + dim.ylength - 1) % dim.ylength) * dim.xlength + (dim.ylength * dim.xlength * z);								 //f(x,y-dy,z)
				N[4][k] = ((x + dim.xlength) % dim.xlength) + ((y + dim.ylength) % dim.ylength) * dim.xlength + (dim.ylength * dim.xlength * ((z + dim.zlength + 1) % dim.zlength)); //f(x,y,z+dz)
				N[5][k] = ((x + dim.xlength) % dim.xlength) + ((y + dim.ylength) % dim.ylength) * dim.xlength + (dim.ylength * dim.xlength * ((z + dim.zlength - 1) % dim.zlength)); //f(x,y,z-dz)
			}
		}
	}

	intvs = 0.25; //interval to print out rho and sum
	intvf = 5;	//interval to print full data

	countt = 0;
	counttf = 0;

	while (t <= finalt)
	{
		sumrho = 0.0;
		sumeta = 0.0;
		S = 0.0;
		SpR = 0.0;
		Sdif = 0.0;

		//parts needed for the bulk free energy term
		for (int a = 0; a < Size; a++)
		{
			etaSqSum = (eta[0][a] * eta[0][a]) + (eta[1][a] * eta[1][a]);
			etaCubeSum = (eta[0][a] * eta[0][a] * eta[0][a]) + (eta[1][a] * eta[1][a] * eta[1][a]);

			rhoCube = rho[a] * rho[a] * rho[a];
			rhoSq = rho[a] * rho[a];
			rhoQuad = rho[a] * rho[a] * rho[a] * rho[a];

			//S = S + rhoCube;

			BFE_partialRho[a] = con.A * ((4 * rhoCube) + ((-(4 * dim.rhoVap) - (4 * dim.rhoSolid) - 2) * rhoSq) + (((4 * dim.rhoSolid * dim.rhoVap) + (2 * dim.rhoVap) + (2 * dim.rhoSolid)) * rho[a]) - (2 * dim.rhoVap * dim.rhoSolid)) + con.B * ((2 * rho[a]) - (6 * rho[a] * etaSqSum) + (4 * rho[a] * etaCubeSum));
			//SpR = SpR +BFE_partialRho[a];

			etaSqZero = eta[0][a] * eta[0][a];
			etaSqOne = eta[1][a] * eta[1][a];
			etaCubeZero = eta[0][a] * eta[0][a] * eta[0][a];
			etaCubeOne = eta[1][a] * eta[1][a] * eta[1][a];

			BFE_partialEta[0][a] = 12 * con.B * (((1 - rho[a]) * eta[0][a]) - ((2 - rho[a]) * etaSqZero) + etaCubeZero + (con.Cgbe * eta[0][a] * etaSqOne / 3));
			BFE_partialEta[1][a] = 12 * con.B * (((1 - rho[a]) * eta[1][a]) - ((2 - rho[a]) * etaSqOne) + etaCubeOne + (con.Cgbe * eta[1][a] * etaSqZero / 3));

			phi = rhoQuad * ((7 * rhoSq) - (18 * rho[a]) + 12);

			Dif[a] = (con.voldif * phi) + (con.vapdif * (1 - phi)) + (con.surfdif * rhoSq * (1 - rho[a]) * (1 - rho[a])) + (con.gbdif * rho[a] * (1 - etaSqSum));
			//Sdif = Sdif +Dif[a];

			//defining h for the rho partial derivative terms

			h[a] = BFE_partialRho[a] - con.surfenergy * (rho[N[0][a]] + rho[N[1][a]] + rho[N[2][a]] + rho[N[3][a]] + rho[N[4][a]] + rho[N[5][a]] - (6 * rho[a]));

			//defining the partial derivative terms with respect to rho and eta

			//for (int a = 0; a<Size; a++)

			//rho(t+1) after the time step the new variable is NewRho. Same goes for eta- eta(t+1)-NewEta

			NewRho[a] = rho[a] + con.deltaT * 0.5 * ((Dif[N[0][a]] * (h[N[0][a]] - h[a])) - (Dif[N[1][a]] * (h[a] - h[N[1][a]])) + (Dif[N[2][a]] * (h[N[2][a]] - h[a])) - (Dif[N[3][a]] * (h[a] - h[N[3][a]])) + (Dif[N[4][a]] * (h[N[4][a]] - h[a])) - (Dif[N[5][a]] * (h[a] - h[N[5][a]])) + (Dif[a] * (h[N[0][a]] + h[N[1][a]] + h[N[2][a]] + h[N[3][a]] + h[N[4][a]] + h[N[5][a]] - (6 * h[a]))));

			NewEta[0][a] = eta[0][a] - con.gbmobility * con.deltaT * (BFE_partialEta[0][a] - (con.gbenergy * (eta[0][N[0][a]] + eta[0][N[1][a]] + eta[0][N[2][a]] + eta[0][N[3][a]] + eta[0][N[4][a]] + eta[0][N[5][a]] - (6 * eta[0][a]))));

			NewEta[1][a] = eta[1][a] - con.gbmobility * con.deltaT * (BFE_partialEta[1][a] - (con.gbenergy * (eta[1][N[0][a]] + eta[1][N[1][a]] + eta[1][N[2][a]] + eta[1][N[3][a]] + eta[1][N[4][a]] + eta[1][N[5][a]] - (6 * eta[1][a]))));
		}

		//restricting field variable values to greater than 0.5

		for (int i = 0; i < Size; i++)
		{ //density field variable
			/* if (rho[i]>1)
        {
        rho[i]=1;
        }*/

			if (rho[i] < 0.0001)
			{

				rho[i] = dim.rhoVap;
			}
			else
			{

				rho[i] = rho[i];
			}
		}

		for (int i = 0; i < Size; i++)
		{ //density field variable
			for (int g = 0; g < dim.noOfparticles; g++)
			{ //orientation parameter
				/*if (eta[g][i]>1)
            {
            eta[g][i]=1;
            }*/
				if (eta[g][i] < 0.0001)
				{
					eta[g][i] = dim.etaVap;
				}
				else
				{
					eta[g][i] = eta[g][i];
				}
			}
		}

		SwitchRho = rho;
		rho = NewRho;
		NewRho = SwitchRho;
		//SwitchRho [a] = rho[a]; rho[a] = NewRho[a]; NewRho[a] = witchRho[a];

		//SwitchEta[0][a] = eta[0][a];
		//SwitchEta[1][a] = eta[1][a];
		for (int g = 0; g < dim.noOfparticles; g++)
		{
			SwitchEta[g] = eta[g];
		}
		//eta[0][a] = NewEta[0][a];
		//eta[1][a] = NewEta[1][a];
		for (int g = 0; g < dim.noOfparticles; g++)
		{
			eta[g] = NewEta[g];
		}
		//NewEta[0][a] = SwitchEta[0][a];
		//NewEta[1][a] = SwitchEta[1][a];
		for (int g = 0; g < dim.noOfparticles; g++)
		{
			NewEta[g] = SwitchEta[g];
		}

		//Enforcing conservation of density

		sumrho = 0.0;

		for (int i = 0; i < Size; i++)
		{ // still has nan
			sumrho = sumrho + rho[i];
		}

		sizebv = 0;
		for (int i = 0; i < Size; i++)
		{
			if (rho[i] > 0.0041)
			{
				sizebv = sizebv + 1;
			}
		}

		intvrho = (sumrhoInit - sumrho) / sizebv;
		for (int i = 0; i < Size; i++)
		{
			if (rho[i] > 0.0041)
			{
				tempr = rho[i];
				rho[i] = tempr + intvrho;
			}
		}
		//cout<<"time = "<<t<<"\t"<<sumrhoInit<<"\t"<<sizebv<<"\t"<<sumrho<<"\n";
		sumrho = 0.0;
		for (int i = 0; i < Size; i++)
		{ // rho still has a nan in it
			sumrho = sumrho + rho[i];
		}

		/// NEWRHO is the problem
		t = t + con.deltaT;


	ofstream fout("Combinations.txt", ios::out | ios::app);
	switch (current)/* This switch statement runs the right combination of parameters*/
	{
		case 1:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0);
				goto R;
			}
			break;

		case 2:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1);
				goto R;
			}
			break;

		case 3:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "," << mapback_dict[i2] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i2], &lowerbound[i2], &upperbound[i2], &MAX[i2], &t);
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1, i2);
				goto R;
			}
			break;
		case 4:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "," << mapback_dict[i2] << "," << mapback_dict[i3] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i3], &lowerbound[i3], &upperbound[i3], &MAX[i3], &t);
			modify(&sumrho, &params[i2], &lowerbound[i2], &upperbound[i2], &MAX[i2], &t);
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1, i2, i3);
				goto R;
			}
			break;

		case 5:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "," << mapback_dict[i2] << "," << mapback_dict[i3] << "," << mapback_dict[i4] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i4], &lowerbound[i4], &upperbound[i4], &MAX[i4], &t);
			modify(&sumrho, &params[i3], &lowerbound[i3], &upperbound[i3], &MAX[i3], &t);
			modify(&sumrho, &params[i2], &lowerbound[i2], &upperbound[i2], &MAX[i2], &t);
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1, i2, i3, i4);
				goto R;
			}
			break;

		case 6:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "," << mapback_dict[i2] << "," << mapback_dict[i3] << "," << mapback_dict[i4] << "," << mapback_dict[i5] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i5], &lowerbound[i5], &upperbound[i5], &MAX[i5], &t);
			modify(&sumrho, &params[i4], &lowerbound[i4], &upperbound[i4], &MAX[i4], &t);
			modify(&sumrho, &params[i3], &lowerbound[i3], &upperbound[i3], &MAX[i3], &t);
			modify(&sumrho, &params[i2], &lowerbound[i2], &upperbound[i2], &MAX[i2], &t);
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1, i2, i3, i4, i5);
				goto R;
			}
			break;

		case 7:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "," << mapback_dict[i2] << "," << mapback_dict[i3] << "," << mapback_dict[i4] << "," << mapback_dict[i5] << "," << mapback_dict[i6] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i6], &lowerbound[i6], &upperbound[i6], &MAX[i6], &t);
			modify(&sumrho, &params[i5], &lowerbound[i5], &upperbound[i5], &MAX[i5], &t);
			modify(&sumrho, &params[i4], &lowerbound[i4], &upperbound[i4], &MAX[i4], &t);
			modify(&sumrho, &params[i3], &lowerbound[i3], &upperbound[i3], &MAX[i3], &t);
			modify(&sumrho, &params[i2], &lowerbound[i2], &upperbound[i2], &MAX[i2], &t);
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1, i2, i3, i4, i5, i6);
				goto R;
			}
			break;
		case 8:
			if (writeonce)
			{
				fout << "Combination: [" << mapback_dict[i0] << "," << mapback_dict[i1] << "," << mapback_dict[i2] << "," << mapback_dict[i3] << "," << mapback_dict[i4] << "," << mapback_dict[i5] << "," << mapback_dict[i6] << "," << mapback_dict[i7] << "]" << endl;
				fout.close();
				writeonce = false;
			}
			modify(&sumrho, &params[i7], &lowerbound[i7], &upperbound[i7], &MAX[i7], &t);
			modify(&sumrho, &params[i6], &lowerbound[i6], &upperbound[i6], &MAX[i6], &t);
			modify(&sumrho, &params[i5], &lowerbound[i5], &upperbound[i5], &MAX[i5], &t);
			modify(&sumrho, &params[i4], &lowerbound[i4], &upperbound[i4], &MAX[i4], &t);
			modify(&sumrho, &params[i3], &lowerbound[i3], &upperbound[i3], &MAX[i3], &t);
			modify(&sumrho, &params[i2], &lowerbound[i2], &upperbound[i2], &MAX[i2], &t);
			modify(&sumrho, &params[i1], &lowerbound[i1], &upperbound[i1], &MAX[i1], &t);
			if (modify(&sumrho, &params[i0], &lowerbound[i0], &upperbound[i0], &MAX[i0], &t, 1))
			{
				reassign(&con, current, i0, i1, i2, i3, i4, i5, i6, i7);
				goto R;
			}
			break;

		} // end switch
	}	 //end while loop

	myfile.close();
	myfile1.close();
	myfile2.close();

	// send files
	//char cmd[] = {"python email_sender.py"};
	//system(cmd);
	return 0;
}
