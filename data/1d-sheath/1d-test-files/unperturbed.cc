#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define IMAX 10000
#define STEP 0.01

int main ()
{
	double Z[IMAX], PHI[IMAX];
	PHI[0] = - 2.838;//775;
	Z[0] = 0.0;

	for (int i=1; i<IMAX; i++)
	{
		Z[i] = i*STEP;
		PHI[i] = PHI[i-1] + STEP * sqrt( 2.0 * sqrt(1.0 - 2.0 * PHI[i-1]) + 2.0 * exp(PHI[i-1]) - 4.0);
	}




	ofstream myfile;
	myfile.open ("unperturbed.dat");
	for (int i=0; i<IMAX; i++)
	{
		myfile << Z[i] << "\t" << PHI[i] << "\t" << 1.0 / sqrt(1.0 - 2.0 * PHI[i]) << "\t" << exp(PHI[i])
		       << "\t" << - sqrt(1.0 - 2.0 * PHI[i]) << "\t"
		       << - sqrt( 2.0 * sqrt(1.0 - 2.0 * PHI[i]) + 2.0 * exp(PHI[i]) - 4.0) << endl;
		// print << position << potential << ion density << electron density << ion velocity << electric field
	}
	myfile.close();

	return 0;
}

