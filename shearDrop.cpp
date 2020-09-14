/*
Main file for simulating an initially circular drop subject to a shear flow in 2D. 

Author: Colton Bryant
*/

//includes
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<string>

//stokes header
#include "stokesGrid.hpp"

//level set headers
#include "LevelSet/level/initialfunc.h"
#include "LevelSet/level/initfuncs/circle.h"
#include "LevelSet/level/initfuncs/linearfunc.h"
#include "fluidInterface.hpp"


using namespace std;
using namespace levelset;

//main function
/**
* Sets up a test problem with an intially circular drop in shear flow. \n
* Inputs: nx, ny - discrete grid size to be used \n
* Output: data written to: \n 
* u_<timestep>.out \n
* v_<timestep>.out \n
* p_<timestep>.out \n
* phi_<timestep>.out
*/
int main(int argc, char* argv[]) {
	
	//user inputs
	int nx 			= atoi(argv[1]); //number of x grid points
	int ny 			= atoi(argv[2]); //number of y grid points

	//set the computational domain
	double xMin = -1.;
	double yMin = -1.;
	double xMax = 1.;
	double yMax = 1.;
	
	//flow variables
	double sigma = 1.;   //surface tension
	double uYMax = yMax; //upper wall velocity
	double uYMin = yMin; //lower wall velocity

	//timesteps
	double tFinal = 1.;
	double dt = 1.e-3;
	int numTSteps = (int) ceil(tFinal/dt);
	int dumpFreq  = 100;  //how often to dump data 	
	
	//data dump base filenames
	string xVelOut = "u";
	string yVelOut = "v";
	string presOut = "p";
	string phiOut  = "phi";

	//create the grids
	stokesGrid background(nx, ny, xMin, xMax, yMin, yMax); //creates the staggered grid object
	fluidInterface drop(nx, ny, xMin, xMax, yMin, yMax, sigma); //creates level set objects
	
	//Stokes flow parameters
	//background.turnOnVerbose(); //option to ask the Stokes solver for residual information as it iterates

	background.setTol(1.e-8);   //tolerance for the Stokes solver 
	background.setUppWallVelocity(uYMax);  //wall velocity at y=ymax
	background.setLowWallVelocity(uYMin); //wall velocity at y=ymin
	background.setOmegaU(0.4); //relaxation parameter for momentum equations
	background.setOmegaP(0.4); //relaxation parameter for pressure update
	
	
	//initialize a circle
	InputParams initParams("circle.in"); //read in parameter for the initial shape (will prompt user if file doesn't exist)
	Circle initCirc(&initParams); //one of many available initial shapes. Others in LevelSet/src/initfuncs/
	initParams.CloseTemplate();	
	

	//set the initial shape
	drop.initialize(initCirc); //initializes the level set function
	
	cout << "Initialization complete, starting temporal evolution...\n"; //let user know where we are
	
	//main temporal loop
	for(int k=0; k<=numTSteps; k++) {
		
		//compute surface tension
		drop.computeSurfaceTension(&background); 
		
		//solve stokes equations
		background.solve();
		
		//check if we should dump data (this is in an odd spot in the loop so flow data is in sync with interface's location)
		if( k%dumpFreq == 0 ) {
			cout << "Writing data at t = " << k*dt << "\n";
			cout << "-----------------------------------------------------------------\n";
			//dump flow data
			background.dumpFlowData(xVelOut+"_"+to_string(k)+".out", yVelOut+"_"+to_string(k)+".out", presOut+"_"+to_string(k)+".out");
			//dump level set data
			drop.dumpSlice("phi", phiOut+"_"+to_string(k)+".out"); //first arg. is name of the slice you want to dump. See fluidInterface.cpp for more
		}
		
		//extend the velocity
		drop.setSpeed(&background);

		//update interface location
		drop.upwindAdvance(dt);
	}
	cout << "Simulation complete. Exiting...\n";
	//success
	return 0;
}


