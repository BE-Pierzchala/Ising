/*
 Program to simulate 2D Ising model

Compiled executable takes arguments in command line in this order:
	- n, the size of nxn lattice
	- t, number of timesteps for the simulation
	- m, number of simulations to average over 
	- random, either 1 = random or 0 = aligned starting grid
	- T1 T2 ..., temperatures   IN UNITS OF J/K_b   separated by space for which to calcule
And returns 4 files, (each temperature separated by \n):
	- Energies.txt, file with energy per site at each time step 
	- EnergiesSTD.txt, standard deviations of each of the above energies
	- Mags.txt, file with magnetisations per site at each time step
	- MagsSTD.txt, standard deviations of the above magnetisations

Short description of algorithm:
	Start with the system with all spins aligned ( random causes problems at low T ).
	At each timestep choose n*n sites at random and decide wheter to flip the spin and
	record new energy and new magnetisation.
	Repeat that m times, at the end obtain averages and STD's at each time step

This program uses functions from myFunctions.h 

example use of executable: ./executable 100 200 50 0 1 2
creates 100 x 100 system, with 200 steps, averaged 50 times, spins orignally aligned
at temperatures 1 and 2
 */


//== Include block ========================================
#include "myFunctions.h"
#include <iostream>
#include <tuple>
#include <fstream>
#include <cmath>
using namespace std;
//=========================================================

int main(int argc, char *argv[]){

//=== Declaring variables =================================    
    int n = stoi( argv[1] );					// get size of the n x n lattice from imput
    int t = stoi( argv[2] );					// get number of 'time steps' 
    double m = stod( argv[3] );					// get number of simulations to run and average
	bool random = stoi( argv [4] );				// Whether grid is to be random or aligned to start with
	double H = stod(argv[5]);					// Get the magnetic field 
	double Temperatures[argc-6];				// Get temperatures
	for( int i = 6; i < argc; i++ ){
		Temperatures[i-6] = stod( argv[i] );
	}
    double E[ (int)m ][t], M[ (int)m ][t];		// Initiate arrays to keep Energies and magnetisation
    int spins[n*n];								// Represent nxn system as 1D n*n array

	// Initialize txt files to which write
	ofstream Energies("Energies.txt"), EnergiesSTD("EnergiesSTD.txt"), Mags("Mags.txt"), MagsSTD("MagsSTD.txt");
//=========================================================    


	for(int k = 0; k < argc - 5; k++){				// Repeat for all temperatures
		double T = Temperatures[k];					// dummy variable to hold current temperature
//=== Simulation ==========================================

		double dE,dM; 								// initialize dummy variables for extracting data from advance(...)
		for( int i = 0; i < m; i++){				// Simulate m times
			
			initialize(spins,n, random);			// Initialize spins grid to be either aligned or random 
			tie(E[i][0],M[i][0]) = getEM(spins, n, random, H);	// get initial values of energy and magnetisation
			for( int j = 1; j < t; j++){			// Update the system t times, at each step recording new values
				tie(dE,dM) = update(spins, n, T,H);	// of energy and magnetisation per site
				E[i][j] = E[i][j-1] + dE/(n*n);
				M[i][j] = M[i][j-1] + dM/(n*n);
			}
			
		}
//=========================================================

//=== Get means, STD and save results =====================

			for( int i = 0; i < t; i++){ 

			// Initalize dummy variables to hold sums of values and squares of values respectively
			double E_sum = 0, EE_sum = 0, M_sum = 0, MM_sum = 0;
		
			for( int j = 0; j < m; j++){			// Sum all simulations at a given time step
				E_sum += E[j][i];					// Sum energies
				EE_sum += E[j][i]*E[j][i];
				
				M_sum +=  M[j][i];					// Sum magnetisations
				MM_sum += M[j][i]*M[j][i];
			}

			Energies<<E_sum/m;										// Write to file mean energy 
			EnergiesSTD<< pow( EE_sum/m - E_sum/m * E_sum/m, 0.5 );	// Write to file STD of mean energy	
			Mags<<M_sum/m;											// -||- magnetisation
			MagsSTD<< pow( MM_sum/m - M_sum/m * M_sum/m, 0.5 );		// -||- magnetisation
			
			if( i != t - 1 ){						// Put ',' between all written elements
				Energies<<",";
				EnergiesSTD<<",";
				Mags<<",";
				MagsSTD<<",";
			}
		}
		if( k != argc - 6){
			Energies<<endl;							// At the end of each temperature put '\n'
			EnergiesSTD<<endl;						// If it's not the last T
			Mags<<endl;
			MagsSTD<<endl;
		}
	}
//=========================================================

return 0;
}
