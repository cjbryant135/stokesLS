/*
Header for fluid/fluid interface class. Used as an interface with Prof. Chopp's level set library.

Author: Colton Bryant
*/

#ifndef FLUIDINTERFACE_H
#define FLUIDINTERFACE_H

//includes
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<map>
#include "LevelSet/level/uniformmesh2d.h"
#include "LevelSet/level/um2boundary.h"
#include "LevelSet/level/um2linear.h"
#include "LevelSet/level/um2xperiodic.h"
#include "LevelSet/level/initialfunc.h"
#include "stokesGrid.hpp"

using namespace std;
using namespace levelset;

/**
* A class to interact with Prof. Chopp's level set library. Used to compute surface tension forces, compute the 
* velocity extension, and update the interface location.
*/
class fluidInterface
{
public:
	//constructor
	fluidInterface(const int nx, const int ny, const double xMin, const double xMax, const double yMin, 
		const double yMax, const double sigma);
	//setters
	void setST(const double sigma);
	void setDeltaWidth(const double eps);
	//getters
	double getXMin();
	double getYMin();
	double getXMax();
	double getYMax();
	double getDeltaWidth();
	//methods
	int nameToIndx(const string slice);
	void setSpeed(const stokesGrid* sGrid);
	void upwindAdvance( const double dt);
	void dumpSlice(const string slice, const string filename);
	void initialize(const InitialFunc& f);
	void computeSurfaceTension();
	
	void computeSurfaceTension(stokesGrid* sGrid); //interpolates surface tension to the staggered grid
	
	double distToInt(levelset::Bicubic* fit, const double xi, const double yj, const double xMin, const double xMax, 
		const double yMin, const double yMax, double &xC, double &yC, bool &clean);
	
	template <typename T> int sgn(T val) {return (T(0) < val) - (val < T(0));}

	~fluidInterface(void); //destructor
private:
	double mSigma, mEps; //surface tension and width of discrete delta used in surface tension
	UniformMesh2D* grid; //level set grid object
	//double slice map
	map<string, int> dIndx;
	//integer slice map
	map<string, int> iIndx;

	int dSlices;
	int iSlices;
	int mNx;
	int mNy;

};
#endif
