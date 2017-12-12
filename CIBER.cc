/*  CIBER.cc  */
/*  Cold Ions with Boltzmann Electron Relation	 */
/*  Some commenting and uncommenting of lines 93-108 is required in order to specify the shape of the sheath surface */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

#define IMAX 50				/* number of cells in x-direction */
#define JMAX 1050 			/* number of cells in y-direction */
#define INV_H 50.0			/* cell width */
#define TOL 1e-6			/* Poisson solver error tolerance */
#define RELAX 1.9			/* relaxation parameter for SOR */
#define TIME_STEP 9999.0		/* time between writing data to file */
#define END_TIME 10.0			/* end time of simulation */
#define PI 3.141592653589793		/* pi */
#define MIME 1836.152663		/* ion-to-electron mass ratio (for grain charging calculation) */

#define FLOW_VELOCITY 1.0		/* (magnitude of) velocity of cold ion flow (set as 1 for Bohm sheath) */
#define GRAIN_POTENTIAL -10.0		/* initial potential of sphere / surface */
#define GRAIN_Z 19.9			/* z-position of grain (bottom of simulation domain is z=0) OR height of sheath */
#define GRAIN_RAD 1.0			/* radius of spherical dust grain */

#define SHEATH				/* turn this on for a sheath simulation (default is flow past objects) */
#define CARTESIAN			/* turn this on for cartesian x-y coordiates (default is axisymmetric r-z) */

/*  Need a grid  */
typedef struct Grid {
	double dt;
	double LS[IMAX+2][JMAX+2];	/* level-set function */
	double E_pot[IMAX+2][JMAX+2];	/* electric potential */
	double ION_U[IMAX+1][JMAX+2];
	double ION_V[IMAX+2][JMAX+1];
	double ION_N[IMAX+2][JMAX+2];
	double TEMP[IMAX+2][JMAX+2];
	double TEMP2[IMAX+2][JMAX+2];
	double GRAIN_PHI;
	double Ii[5];
	double DRAG[5];
} Grid;


/*  Define a function to initialize the grid values  */
void init_grid( Grid* grid ) {
	grid->GRAIN_PHI = GRAIN_POTENTIAL;
	for(int surf=0; surf<5; surf++) {
		grid->Ii[surf] = 0.0;
		grid->DRAG[surf] = 0.0;
	}
	double R, Z, UN, HALF;
	UN = sqrt( FLOW_VELOCITY*FLOW_VELOCITY - 2.0*grid->GRAIN_PHI);
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
#ifdef SHEATH
			Z = GRAIN_Z;
			grid->E_pot[i][j] = -6.0/(Z*Z);
			grid->ION_U[i][j] = 0.0;
			grid->ION_V[i][j] = -FLOW_VELOCITY - 6.0/(Z*Z);
#else
			grid->E_pot[i][j] = 0.0;
			/* Initialize velocity so that upwinding points away from grain - use solution for potential flow with 
			   sink (normal velocity -1) at radius 1 */
			R = (double(i)-0.5)/INV_H;
			Z = (double(j)-GRAIN_Z*INV_H-0.5)/INV_H;
			HALF = 0.5/INV_H;
			if ( i!=IMAX+1) {
				grid->ION_U[i][j] = -(R+HALF)*UN*GRAIN_RAD*GRAIN_RAD/((R+HALF)*(R+HALF)+Z*Z)
											/sqrt((R+HALF)*(R+HALF)+Z*Z)
						+ 1.5*GRAIN_RAD*GRAIN_RAD*GRAIN_RAD*Z*(R+HALF)*FLOW_VELOCITY
						/((R+HALF)*(R+HALF)+Z*Z)/((R+HALF)*(R+HALF)+Z*Z)/sqrt((R+HALF)*(R+HALF)+Z*Z);
			}
			if ( j!=JMAX+1) {
				grid->ION_V[i][j] = -(Z+HALF)*UN*GRAIN_RAD*GRAIN_RAD/(R*R+(Z+HALF)*(Z+HALF))
											/sqrt(R*R+(Z+HALF)*(Z+HALF))
						- (1.0 + 0.5*(R*R-2.0*(Z+HALF)*(Z+HALF))*GRAIN_RAD*GRAIN_RAD*GRAIN_RAD
							/(R*R+(Z+HALF)*(Z+HALF))/(R*R+(Z+HALF)*(Z+HALF))
										/sqrt(R*R+(Z+HALF)*(Z+HALF)))*FLOW_VELOCITY;
			}
#endif
			grid->TEMP[i][j] = grid->ION_U[i][j];
			grid->TEMP2[i][j] = grid->ION_V[i][j];
		}
	}

	// Set up LS function
	double r,z;
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
#ifdef SHEATH
			// Sinusoidal electrode
//			grid->LS[i][j] = -( (double(j)-0.5)/INV_H - 1.0 )
//						 + 0.1*cos(3.14159*(double(i)-0.5)/INV_H);

			// square trench
			if (i<IMAX/2) {
				grid->LS[i][j] = 0.1 - (double(j)-0.5)/INV_H;
			} else {
				grid->LS[i][j] = 1.1 - (double(j)-0.5)/INV_H;
			}

			// half-ellipsoid
//			double half_width = 0.2;
//			double half_height = 0.0;
//			r = ((double(i)-0.5)/INV_H);
//			z = ((double(j)-0.5)/INV_H - 0.1);
//			grid->LS[i][j] = max( 0.1 - sqrt(r*r/half_width/half_width + z*z/half_height/half_height), -z );
#else
			// Sphere
			r = ((double(i)-0.5)/INV_H);
			z = ((double(j)-0.5)/INV_H - GRAIN_Z);
			grid->LS[i][j] = GRAIN_RAD - sqrt(r*r+z*z);
#endif
		}
	}

	// Ion velocity/density zero inside object (so no density flux leaves it)
	for (int i=0; i<=IMAX+1; i++) {
		for (int j=0; j<=JMAX+1; j++) {
			grid->ION_N[i][j] = 0.0;
#ifdef SHEATH
			Z = GRAIN_Z;
			if (grid->LS[i][j] < 0.0) grid->ION_N[i][j] = 1.0 - 6.0/(Z*Z);
#else
			if (grid->LS[i][j] < 0.0) grid->ION_N[i][j] = 1.0;	
#endif
		}
	}
}


/*  Define a function to write grid-centred values to file labelled with time t */
void write_grid( Grid* grid, double t) {
	ostringstream fileNameStream("");
#ifdef SHEATH
	fileNameStream << "output/t" << t << "_gridvalues.dat";
#else
	fileNameStream << "output/U" << FLOW_VELOCITY << "A" << GRAIN_RAD << "t" << t << "_gridvalues.dat";
#endif
	string fileName = fileNameStream.str();
	ofstream myfile;
	myfile.open ( fileName.c_str() );

	double Er, Ez;

	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			/* Reading out x / y / level-set fn / electric potential / Er / Ez / U / V / N / Energy */
			Er = -(grid->E_pot[i+1][j]-grid->E_pot[i-1][j])*0.5*INV_H;
			Ez = -(grid->E_pot[i][j+1]-grid->E_pot[i][j-1])*0.5*INV_H;
			myfile << (float(i)-0.5)/INV_H << "\t" << (float(j)-0.5)/INV_H << "\t" 
				<< grid->LS[i][j] << "\t" << grid->E_pot[i][j] << "\t"
				<< Er << "\t" << Ez << "\t" << 0.5*(grid->ION_U[i][j]+grid->ION_U[i-1][j]) << "\t"
				<< 0.5*(grid->ION_V[i][j]+grid->ION_V[i][j-1]) << "\t" << grid->ION_N[i][j] << "\t"
				<< 0.5*(0.25*(grid->ION_V[i][j]+grid->ION_V[i][j-1])*(grid->ION_V[i][j]+grid->ION_V[i][j-1])
					+ 0.25*(grid->ION_U[i][j]+grid->ION_U[i-1][j])*(grid->ION_U[i][j]+grid->ION_U[i-1][j]) )
								+ grid->E_pot[i][j] << endl;
		}
	myfile << endl;
	}
	myfile.close();
}

/*  And a function to write timeseries of energy errors (average and maximum) */
void write_energy( Grid* grid, double t) {
	double energy=0.0, maxenergy=0.0, xmax=0.0, ymax=0.0, avenergy=0.0, npoints=0.0;

	ostringstream fileNameStream("");
#ifdef SHEATH
	fileNameStream << "output/energy.dat";
#else
	fileNameStream << "output/U" << FLOW_VELOCITY << "A" << GRAIN_RAD << "energy.dat";
#endif
	string fileName = fileNameStream.str();
	ofstream myfile2;
	myfile2.open ( fileName.c_str(), ios_base::app);

	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] <  0.0) {
				energy = 0.5*(0.25*(grid->ION_V[i][j]+grid->ION_V[i][j-1])*(grid->ION_V[i][j]+grid->ION_V[i][j-1])
					+ 0.25*(grid->ION_U[i][j]+grid->ION_U[i-1][j])*(grid->ION_U[i][j]+grid->ION_U[i-1][j]) )
								+ grid->E_pot[i][j];
				energy = abs(energy-0.5);
				if (energy > maxenergy) {
					maxenergy = energy;
					xmax = (double(i)-0.5)/INV_H;
					ymax = (double(j)-0.5)/INV_H;
				}
				avenergy += energy;
				npoints += 1.0;
			}
		}
	}
	avenergy /= npoints;

	myfile2 << t << "\t" << avenergy << "\t" << maxenergy << "\t" << xmax << "\t" << ymax << endl;
	myfile2.close();
}


/*  And a function to write timeseries of charge and drag force */
void write_phi_and_drag( Grid* grid, double t) {
	ostringstream fileNameStream("");
	fileNameStream << "output/macroscopicU" << FLOW_VELOCITY << "A" << GRAIN_RAD << ".dat";
	string fileName = fileNameStream.str();
	ofstream myfile3;
	myfile3.open ( fileName.c_str(), ios_base::app);
	myfile3 << t*TIME_STEP << "\t" << grid->GRAIN_PHI
				<< "\t" << grid->Ii[0] << "\t" << grid->Ii[1] << "\t" << grid->Ii[2] 
				<< "\t" << grid->Ii[3] << "\t" << grid->Ii[4] 
				<< "\t" << grid->DRAG[0] << "\t" << grid->DRAG[1] << "\t" << grid->DRAG[2]
				<< "\t" << grid->DRAG[3] << "\t" << grid->DRAG[4] << endl;
	myfile3.close();
}


/*  Function to solve Poisson-Boltzmann equation in plasma region */
void solve_PB_eqn( Grid* grid) {
	double temp, err, maxerr, R;			/* Use min norm for residual and tolerance */
	double S, E, W;					/* Switches to zero for points on the boundary */

#ifndef SHEATH
	for (int i=0; i<=IMAX+1; i++) {
		grid->E_pot[i][JMAX+1] = 0.0;
	}
#endif
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] >= 0.0) {
				grid->E_pot[i][j] = grid->GRAIN_PHI;
			}
		}
	}
	do {
		// Successive-over-relaxation.
		maxerr = 0.0;
		for (int i=1; i<=IMAX; i++) {
			for (int j=1; j<=JMAX; j++) {
				if (grid->LS[i][j] <  0.0) {	// only solve in vacuum
					E = 1.0;
					W = 1.0;
					S = 1.0;
					if (i == 1) W = 0.0;
					if (i == IMAX) E = 0.0;
					if (j == 1) S = 0.0;

#ifdef CARTESIAN
					temp = (1.0-RELAX)*grid->E_pot[i][j] + RELAX *
						( E*grid->E_pot[i+1][j] + W*grid->E_pot[i-1][j]
							+ grid->E_pot[i][j+1] + S*grid->E_pot[i][j-1] 
							+ (grid->ION_N[i][j] - exp(grid->E_pot[i][j]))/INV_H/INV_H )
							 	/(1.0+S+E+W) ;
#else
					R = double(i)-0.5;
					temp = (1.0-RELAX)*grid->E_pot[i][j] + RELAX *
						( (R+0.5)*E*grid->E_pot[i+1][j] + (R-0.5)*W*grid->E_pot[i-1][j]
							+ R*grid->E_pot[i][j+1] + R*S*grid->E_pot[i][j-1] 
							+ R*(grid->ION_N[i][j] - exp(grid->E_pot[i][j]))/INV_H/INV_H )
							 	/((1.0+S)*R+(R+0.5)*E+(R-0.5)*W) ;
#endif
					err = abs(grid->E_pot[i][j] - temp);
					if (err > maxerr) maxerr = err;
					grid->E_pot[i][j] = temp;
				}
			}
		}
	} while (maxerr > TOL);

	for (int j=1; j<=JMAX-1; j++) {
		grid->E_pot[0][j] = grid->E_pot[1][j];
		grid->E_pot[IMAX+1][j] = grid->E_pot[IMAX][j];
	}
#ifndef SHEATH
	for (int i=0; i<=IMAX+1; i++) {
		grid->E_pot[i][0] = grid->E_pot[i][1];
	}
#endif
}


/*  This function sets DT to give numerically stable calculations according to CFL condition */
void set_dt( Grid* grid) {
	double max_velocity = 0.0;
	for (int i=1; i<=IMAX; i++) {			// loop such that all interior u and v covered
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] <  0.0) {
				max_velocity = max( abs(grid->ION_U[i][j]), max_velocity);
				max_velocity = max( abs(grid->ION_V[i][j]), max_velocity);
			}
		}
	}
	grid->dt = 0.25 /(INV_H*max_velocity);
}


/*  Function to calculate ion and electron currents to spherical grain and update GRAIN_PHI */
bool charge_grain( Grid* grid) {
	/* One-way Maxwellian flux of electrons */
#ifdef CARTESIAN
	double Ie = - 2.0 * PI * GRAIN_RAD * sqrt(0.5*MIME/PI) * exp(grid->GRAIN_PHI);
#else
	double Ie = - 4.0 * PI * GRAIN_RAD * GRAIN_RAD * sqrt(0.5*MIME/PI) * exp(grid->GRAIN_PHI);
#endif
	double Ii, DRAG, R, Z, N, U, V, Phi, Fr, Fz;
	int I, J1, J2;

	for(int surf=0; surf<5; surf++) {
		/* Current and drag force integrals on cylindrical surfaces (Set to 5 different surfaces) */
		if (surf==0) {
			I = int( (GRAIN_RAD + 0.5) * INV_H );
			J1 = int( (GRAIN_Z - GRAIN_RAD - 0.5) * INV_H );
			J2 = int( (GRAIN_Z + GRAIN_RAD + 0.5) * INV_H );
		} else if (surf==4) {
			I = IMAX-2;		// FIXME: place at 0.5 from edge??
			J1 = 2;
			J2 = JMAX-2;
		} else {
			I = surf*IMAX/4 + (4-surf)*int( GRAIN_RAD * INV_H )/4;
			J1 = (4-surf)*int( (GRAIN_Z - GRAIN_RAD) * INV_H )/4;
			J2 = surf*JMAX/4 + (4-surf)*int( (GRAIN_Z + GRAIN_RAD) * INV_H )/4;
		}
		Ii = 0.0;
		DRAG = 0.0;
		for (int i = 1; i <= I; i++) {
			N = 0.5*(grid->ION_N[i][J1] + grid->ION_N[i][J1+1]);
			V = grid->ION_V[i][J1];
			Phi = 0.5*(grid->E_pot[i][J1] + grid->E_pot[i][J1+1]);
			Fr = 0.25*(grid->E_pot[i+1][J1] + grid->E_pot[i+1][J1+1]
						- grid->E_pot[i-1][J1] - grid->E_pot[i-1][J1+1])*INV_H;
			Fz = (grid->E_pot[i][J1+1] - grid->E_pot[i][J1])*INV_H;
#ifdef CARTESIAN
			Ii += N * V / INV_H;
			DRAG += ( 0.5*Fz*Fz - 0.5*Fr*Fr - N*V*V - exp(Phi) ) / INV_H;
#else
			R = (double(i)-0.5) / INV_H;
			Ii += 2.0 * PI * R * N * V / INV_H;
			DRAG += 2.0 * PI * R * ( 0.5*Fz*Fz - 0.5*Fr*Fr - N*V*V - exp(Phi) ) / INV_H;
#endif

			N = 0.5*(grid->ION_N[i][J2] + grid->ION_N[i][J2+1]);
			V = grid->ION_V[i][J2];
			Phi = 0.5*(grid->E_pot[i][J2] + grid->E_pot[i][J2+1]);
			Fr = 0.5*(grid->E_pot[i+1][J2] + grid->E_pot[i+1][J2+1]
						- grid->E_pot[i-1][J2] - grid->E_pot[i-1][J2+1])*INV_H;
			Fz = (grid->E_pot[i][J2+1] - grid->E_pot[i][J2])*INV_H;
#ifdef CARTESIAN
			Ii -= N * V / INV_H;
			DRAG -= ( 0.5*Fz*Fz - 0.5*Fr*Fr - N*V*V - exp(Phi) ) / INV_H;
#else
			R = (double(i)-0.5) / INV_H;
			Ii -= 2.0 * PI * R * N * V / INV_H;
			DRAG -= 2.0 * PI * R * ( 0.5*Fz*Fz - 0.5*Fr*Fr - N*V*V - exp(Phi) ) / INV_H;
#endif
		}
		for (int j = J1+1; j <= J2; j++) {
			N = 0.5*(grid->ION_N[I][j] + grid->ION_N[I+1][j]);
			U = grid->ION_U[I][j];
			V = 0.25*(grid->ION_V[I][j] + grid->ION_V[I+1][j] + grid->ION_V[I][j-1] + grid->ION_V[I+1][j-1]);
			Fr = (grid->E_pot[I+1][j] - grid->E_pot[I][j])*INV_H;
			Fz = 0.25*(grid->E_pot[I][j+1] + grid->E_pot[I+1][j+1]
						- grid->E_pot[I][j-1] - grid->E_pot[I+1][j-1])*INV_H;
#ifdef CARTESIAN
			Ii -= N * U / INV_H;
			DRAG -= (Fr * Fz - N * U * V) / INV_H;
#else
			R = double(I) / INV_H;
			Ii -= 2.0 * PI * R * N * U / INV_H;
			DRAG -= 2.0 * PI * R * (Fr * Fz - N * U * V) / INV_H;
#endif
		}
		grid->Ii[surf] = Ii;
		grid->DRAG[surf] = DRAG;
	}

	/* Use inner current to update grain potential - note that vacuum capactance is assumed for sphere...
           And for a cylinder this is absolute garbage, corresponding to ln(b/a) = 1/(2R) */
	grid->GRAIN_PHI += grid->dt*( Ie + grid->Ii[2] )*0.25/(PI * GRAIN_RAD);
}


/*  Function to calculate ion density with second order upwind differences */
bool calc_ion_n( Grid* grid) {

	double R, DNURDR, DNVDZ, RFLUXPLUS, RFLUXMINUS, ZFLUXPLUS, ZFLUXMINUS;
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			if (grid->LS[i][j] <  0.0) {	// only solve in vacuum
				RFLUXPLUS = max(grid->ION_U[i][j],0.0)*grid->ION_N[i][j]
								+ min(grid->ION_U[i][j],0.0)*grid->ION_N[i+1][j];
				RFLUXMINUS = max(grid->ION_U[i-1][j],0.0)*grid->ION_N[i-1][j]
								+ min(grid->ION_U[i-1][j],0.0)*grid->ION_N[i][j];

				ZFLUXPLUS = max(grid->ION_V[i][j],0.0)*grid->ION_N[i][j]
								+ min(grid->ION_V[i][j],0.0)*grid->ION_N[i][j+1];
				ZFLUXMINUS = max(grid->ION_V[i][j-1],0.0)*grid->ION_N[i][j-1]
								+ min(grid->ION_V[i][j-1],0.0)*grid->ION_N[i][j];

				DNVDZ = ZFLUXPLUS - ZFLUXMINUS;

#ifdef CARTESIAN
				DNURDR = RFLUXPLUS - RFLUXMINUS;
				grid->TEMP[i][j] = grid->ION_N[i][j] - grid->dt*INV_H*(DNURDR + DNVDZ);
#else
				R = double(i)-0.5;
				DNURDR = RFLUXPLUS*(R+0.5) - RFLUXMINUS*(R-0.5);
				grid->TEMP[i][j] = grid->ION_N[i][j] - grid->dt*INV_H*(DNURDR/R + DNVDZ);
#endif
			}
		}
	}

	/* FIXME: Is this section really necessary??? */
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX; j++) {
			grid->ION_N[i][j] = min(grid->TEMP[i][j], 1000.0);
			if (grid->LS[i][j] >= 0.0) grid->ION_N[i][j] = 0.0;
			if ( grid->ION_N[i][j] < 0.0 ) grid->ION_N[i][j] = 0.25*(grid->ION_N[i+1][j]+grid->ION_N[i-1][j]
											+grid->ION_N[i][j+1]+grid->ION_N[i][j-1]);
			if ( isnan(grid->ION_N[i][j]) == true ) {
				cout << "N is nan at (" << (double(i)-0.5)/INV_H << ", " << (double(j)-0.5)/INV_H << ")" << endl; 
				return true;
			}
		}
	}

	/* tidy up boundaries (fixed inflow, quasineutral Boltzmann approximation at outer/bottom edges)*/
#ifndef SHEATH
	for (int i=1; i<=IMAX; i++) {
		grid->ION_N[i][0] = exp(grid->E_pot[i][0]);
		grid->ION_N[i][JMAX+1] = 1.0;
	}
#endif
	for (int j=0; j<=JMAX+1; j++) {
		grid->ION_N[0][j] = grid->ION_N[1][j];
#if defined(CARTESIAN) || defined(SHEATH)
		grid->ION_N[IMAX+1][j] = grid->ION_N[IMAX][j];
#else
		grid->ION_N[IMAX+1][j] = exp(grid->E_pot[IMAX+1][j]);
#endif
	}

	return false;
}


/*  Function to calculate ion velocity (convective derivative calculation taken from iEHD) */
bool calc_ion_uv( Grid* grid) {
	double UDUDR, VDUDZ, UDVDR, VDVDZ, DPHIDR, DPHIDZ, DUMINUS, DUPLUS, DVMINUS, DVPLUS;
	double max_vel = sqrt(FLOW_VELOCITY*FLOW_VELOCITY-2.0*grid->GRAIN_PHI);
	for (int i=1; i<=IMAX-1; i++) {				/* Looping over U points only */
		for (int j=1; j<=JMAX; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i+1][j] < 0.0 ) {	// only solve in vacuum

				DUMINUS = (3.0*grid->ION_U[i][j] - 4.0*grid->ION_U[i-1][j] + grid->ION_U[i-2][j])*0.5*INV_H;
				if ( i==1 ) DUMINUS = (grid->ION_U[1][j]-grid->ION_U[0][j])*INV_H;
				DUPLUS = (-grid->ION_U[i+2][j] + 4.0*grid->ION_U[i+1][j] - 3.0*grid->ION_U[i][j])*0.5*INV_H;
				if ( i==IMAX-1 ) DUPLUS = (grid->ION_U[IMAX][j]-grid->ION_U[IMAX-1][j])*INV_H;
				UDUDR = max(grid->ION_U[i][j],0.0)*DUMINUS + min(grid->ION_U[i][j],0.0)*DUPLUS;

				DUMINUS = (3.0*grid->ION_U[i][j] - 4.0*grid->ION_U[i][j-1] + grid->ION_U[i][j-2])*0.5*INV_H;
				if ( j==1 ) DUMINUS = (grid->ION_U[i][1]-grid->ION_U[i][0])*INV_H;
				DUPLUS = (-grid->ION_U[i][j+2] + 4.0*grid->ION_U[i][j+1] - 3.0*grid->ION_U[i][j])*0.5*INV_H;
				if ( j==JMAX ) DUPLUS = (grid->ION_U[i][JMAX+1]-grid->ION_U[i][JMAX])*INV_H;
				VDUDZ = max((grid->ION_V[i][j]+grid->ION_V[i+1][j]+grid->ION_V[i][j-1]
								+grid->ION_V[i+1][j-1])*0.25,0.0) * DUMINUS
					+ min((grid->ION_V[i][j]+grid->ION_V[i+1][j]+grid->ION_V[i][j-1]
								+grid->ION_V[i+1][j-1])*0.25,0.0) * DUPLUS;

				DPHIDR = ( grid->E_pot[i+1][j] - grid->E_pot[i][j] ) * INV_H;

				grid->TEMP[i][j] = grid->ION_U[i][j] - grid->dt*( UDUDR + VDUDZ + DPHIDR );
			}
		}
	}

	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i][j+1] < 0.0 ) {	// only solve in vacuum

				DVMINUS = (3.0*grid->ION_V[i][j] - 4.0*grid->ION_V[i][j-1] + grid->ION_V[i][j-2])*0.5*INV_H;
				if ( j==1 ) DVMINUS = (grid->ION_V[i][1]-grid->ION_V[i][0])*INV_H;
				DVPLUS = (-grid->ION_V[i][j+2] + 4.0*grid->ION_V[i][j+1] - 3.0*grid->ION_V[i][j])*0.5*INV_H;
				if ( j==JMAX-1 ) DVPLUS = (grid->ION_V[i][JMAX]-grid->ION_V[i][JMAX-1])*INV_H;
				VDVDZ = max(grid->ION_V[i][j],0.0)*DVMINUS + min(grid->ION_V[i][j],0.0)*DVPLUS;

				DVMINUS = (3.0*grid->ION_V[i][j] - 4.0*grid->ION_V[i-1][j] + grid->ION_V[i-2][j])*0.5*INV_H;
				if ( i==1 ) DVMINUS = (grid->ION_V[1][j]-grid->ION_V[0][j])*INV_H;
				DVPLUS = (-grid->ION_V[i+2][j] + 4.0*grid->ION_V[i+1][j] - 3.0*grid->ION_V[i][j])*0.5*INV_H;
				if ( i==IMAX ) DVPLUS = (grid->ION_V[IMAX+1][j]-grid->ION_V[IMAX][j])*INV_H;
				UDVDR = max((grid->ION_U[i][j]+grid->ION_U[i][j+1]+grid->ION_U[i-1][j+1]
								+grid->ION_U[i-1][j])*0.25,0.0) * DVMINUS
					+ min((grid->ION_U[i][j]+grid->ION_U[i][j+1]+grid->ION_U[i-1][j+1]
								+grid->ION_U[i-1][j])*0.25,0.0) * DVPLUS;
	
				DPHIDZ = ( grid->E_pot[i][j+1] - grid->E_pot[i][j] ) * INV_H;

				grid->TEMP2[i][j] = grid->ION_V[i][j] - grid->dt*( VDVDZ + UDVDR + DPHIDZ );
			}
		}
	}

	// FIXME: is this bit really necessary?
	for (int i=1; i<=IMAX-1; i++) {				/* Looping over U points only */
		for (int j=1; j<=JMAX; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i+1][j] < 0.0 ) {
				grid->ION_U[i][j] = grid->TEMP[i][j];
				if ( isnan(grid->ION_U[i][j]) == true ) {
					cout << "U is nan at (" << (double(i)-0.5)/INV_H << ", "
									<< (double(j)-0.5)/INV_H << ")" << endl; 
					return true;
				}
			}
		}
	}

	// FIXME: or this?
	for (int i=1; i<=IMAX; i++) {
		for (int j=1; j<=JMAX-1; j++) {
			if ( grid->LS[i][j] <  0.0 || grid->LS[i][j+1] < 0.0 ) {
				grid->ION_V[i][j] = grid->TEMP2[i][j];
				if ( isnan(grid->ION_V[i][j]) == true ) {
					cout << "V is nan at (" << (double(i)-0.5)/INV_H << ", "
									<< (double(j)-0.5)/INV_H << ")" << endl; 
					return true;
				}
			}
		}
	}

	/* top/bottom boundary conditions */
#ifndef SHEATH
	for (int i=0; i<=IMAX; i++) {
		grid->ION_U[i][JMAX+1] = 0.0;
		grid->ION_U[i][0] = grid->ION_U[i][1];
		grid->ION_V[i][0] = -FLOW_VELOCITY + grid->E_pot[i][1]/FLOW_VELOCITY;
		grid->ION_V[i][JMAX] = -FLOW_VELOCITY;
	}
	grid->ION_V[IMAX+1][0] = -FLOW_VELOCITY + grid->E_pot[IMAX+1][1]/FLOW_VELOCITY;
#endif
#if defined(CARTESIAN) || defined(SHEATH)
	for (int j=1; j<=JMAX; j++) {
		grid->ION_U[0][j] = 0.0;
		grid->ION_U[IMAX][j] = 0.0;
		grid->ION_V[0][j] = grid->ION_V[1][j];
		grid->ION_V[IMAX+1][j] = grid->ION_V[IMAX][j];
	}
#else
	double V;
	for (int j=1; j<=JMAX; j++) {
		grid->ION_U[0][j] = 0.0;
		V = 0.5*( grid->ION_V[IMAX][j-1] + grid->ION_V[IMAX][j] );
		grid->ION_U[IMAX][j] = (double(IMAX)-1.0)*grid->ION_U[IMAX-1][j]/double(IMAX)
					+ (double(IMAX)-0.5)*(1.0/V-V )
							*0.5*(grid->E_pot[IMAX][j+1] - grid->E_pot[IMAX][j-1])/double(IMAX);
		grid->ION_V[0][j] = grid->ION_V[1][j];
		grid->ION_V[IMAX+1][j] = -FLOW_VELOCITY + 0.5*(grid->E_pot[IMAX+1][j] + grid->E_pot[IMAX+1][j+1])/FLOW_VELOCITY;
	}
#endif
	
	return false;
}


int main() {
	// set up grid at t=0
	double t = 0.0;
	double integer_t = 0.0;
	double energy_t = 0.0;
	Grid grid;
	init_grid(&grid);
	// get electric potential set up nicely before first file write
	solve_PB_eqn(&grid);
	// Get ready to run...
	bool nan_detected = false;

	// main loop
	do {
		set_dt(&grid);
		if ( t >= integer_t* TIME_STEP ) {
			write_grid(&grid, integer_t);	// NB real time is t, integer_t will have rounding error
#ifndef SHEATH
			write_phi_and_drag(&grid, integer_t);
#endif
			integer_t += 1.0;		// but integer_t is better for filenames
		}
#ifdef SHEATH
		if ( t >= energy_t ) {
			write_energy(&grid, t);
			energy_t += 0.1;
		}
#endif

		// order so that scheme is fully explicit, ie. always use previous grid values to predict next step
		solve_PB_eqn(&grid);
#ifndef SHEATH
		charge_grain(&grid);
#endif
		nan_detected = calc_ion_n(&grid);
		nan_detected = max( calc_ion_uv(&grid), nan_detected);
		t += grid.dt;
		cout << t << "\t" << grid.dt << endl;
	} while (nan_detected == false && t < END_TIME + grid.dt);

	write_grid(&grid, integer_t);
	return 0;
}
