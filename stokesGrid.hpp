/*
Solves Stokes equations on a MAC grid subject to a given body force term. 


Author: Colton Bryant 
*/
#ifndef STOKESGRID_H
#define STOKESGRID_H

//includes
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

using namespace std;

class stokesGrid
{
static const int MAXITS = 1e5; //max iterations for the stokes solver
public:
	//constructor
	stokesGrid(const int nx, const int ny, const double xMin, const double xMax, const double yMin, 
		const double yMax);	
	
	//methods
	//wall velocity setters
	void setUppWallVelocity(const double val);
	void setLowWallVelocity(const double val);

	//setters for relaxation parameters
	void setOmegaU(const double val); //set relaxation parameter for momentum equation
	void setOmegaP(const double val); //set relaxation parameter for pressure equation
	
	void setTol(const double val); //sets tolerance for the solver

	//functions to locate grid points
	double xC(const int i) const;
	double yC(const int j) const;
	double xF(const int i) const;
	double yF(const int j) const;
	
	//flow data access functions
	double u(const int i, const int j) const;
	double v(const int i, const int j) const;
	double p(const int i, const int j) const;
	
	int uI(const int i, const int j) const;
	int vI(const int i, const int j) const;
	int pI(const int i, const int j) const;

	//body force access
	double fx(const int i, const int j);
	double fy(const int i, const int j);

	//body force setters
	void setFX(const int i, const int j, const double val);
	void setFY(const int i, const int j, const double val);

	//solver
	void solve();
	void makeResidualsHuge();
	void zeroOutResiduals();
	bool isConverged();

	//interpolation routines
	double uB(const double x, const double y) const;
	double vB(const double x, const double y) const;
	double pB(const double x, const double y) const;
	double bilin(const double a11, const double a12, const double a21, const double a22, const double ifrac, const double jfrac) const;


	//data dump
	void dumpFlowData(string uOut, string vOut, string pOut);
	void dumpFX(string fxOut);
	void dumpFY(string fyOut);
	
	//turn on verbose output
	void turnOnVerbose();
	//destructor
	~stokesGrid(void);
private:
	double* flowData; //array for velocity and pressure data
	double* bodyXData; //x component of body force terms
	double* bodyYData; //y component of body force terms
	
	//grid paramters
	double mXMin, mXMax, mYMin, mYMax;
	double dx, dy;
	//grid size
	int mNx;
	int mNy;
	
	double tol; //tolerance for the solver
	//relaxation parameters
	double omegaU, omegaP;
	
	double maxURes, maxVRes, maxPRes; //residuals
	double uLow, uHigh; //wall velocities at y=ymin and y=ymax respectively

	int uStart;
	int vStart;
	int pStart;

	bool verbose;
};
#endif
