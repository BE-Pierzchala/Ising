/*
 A header file with functions used in main.cpp of the Ising model project
 Contains functions:
	initialize
	get_dE
	getEM
	update
*/
#ifndef myFunctions_h
#define myFunctions_h

//=== Include block =============
#include <tuple>
using namespace std;
//===============================

void initialize(int spins[], int n, bool random){
	
	// Function to initialize spins array, if random == true the spins will be randomized, otherwise they are all 1 

	if( random ){
		for( int i = 0; i < n*n; i++){
        	spins[i]  =  2*(rand()%2) - 1; 				
		}
	}else{
		for( int i = 0; i < n*n; i++ ){
			spins[i] = 1;
		}
	}
}

double get_dE(int spins [],int position, int n, double H){

	// Function to calculate and return energy difference if a spin at a position in a n x n spins lattice flips ( up -> down or down -> up )    
	// Result is returned in units of J, the interaction energy

	// Define positions of the closet neighbours on the 2D lattice written in a 1D array
    int left,right,up,down;
   
	// Check for all neighbours that they are in the the array range
	// For left and right additionally ensure that the neighbour is in the same row
 
    left = position - 1;
    if( left < 0 || left/n != position/n){
        left += n;
    }
    right = position + 1; 
	if( right > n*n || right/n != position/n){
        right -= n;
    }
    up = position - n;
    if( up < 0 ){
         up +=n*n;
    } 
    down = position + n;
    if( down >= n*n ){
        down -= n*n;
    }
	
	// Return the change in energy   
    return 2*(spins[left] + spins[up] + spins[right] + spins[down])*spins[position] + 2*H*spins[position];
}

tuple<double,double> getEM( int spins[], int n, bool random,double H){

	// Function to get energy and magnetisation of the initial array of spins
	
	if( random ){							// if it's random summ respective parameters
		double E = 0, M = 0;

		for( int i = 0; i < n*n; i++){
			M += spins[i];
			E += -get_dE(spins, i, n,H)/4;
		}
		return make_tuple(E/(n*n),M/(n*n));
	}else{									// if it's aligned return known values 
		return make_tuple(-2-H,1);
	}
}

tuple <double,double> update(int spins[],int n, double T, double H){
    
	/* Function that advances given 1D spins array of size nxn at temperature T and returns total change in energy and magnetisation as a tuple (dE, dM)
	   This approach chooses randomly n*n sites to update, update:
			if change in interaction energy for flipping spin is negative, always flip
			if change in interaction energy for flipping spin is positive, flip only if exp( -dE/(k_b T) ) > random number <1,0>
			else don't flip the spin
		Temperature is provided in units of J/k_b where J is interaction energy coupling
	*/
	
	// Initialize counters of total change in energy and magnetisation
	int dE_tot = 0;
    int dMag = 0;

    for( int i = 0; i < n*n; i++){
        
		double r =  (double)rand()/(double)RAND_MAX; 		// Get a random number in <0,1>
        int pos = (int) (r*n*n);							// Choose a random position to update
        r =  (double)rand()/(double)RAND_MAX; 				// For some reason r has to be drawn again
        double dE = get_dE(spins,pos,n,H);					// get change of interaction energy for flipping spin at pos
        
		// Update, if a flip occurs increment change in energy and magnetisation
        if( dE < 0 ){
            
			spins[pos] *= -1;
            dE_tot += dE;
            dMag += 2*spins[pos];

        }else if( exp( -dE/(double)T ) > double(r)){
            
			spins[pos] *= -1;
            dE_tot += dE;
            dMag += 2*spins[pos];

        }
    }
	
	// return total changes as a tuple
    return make_tuple(dE_tot, dMag);
}


#endif 
