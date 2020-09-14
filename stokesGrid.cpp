/*
Definitions for the stokes grid class.

Notes: It is assumed that boundary conditions are periodic in x and uses no slip/no penetration conditions in y. Non-zero wall velocites can be
specified at y=ymin and y=ymax as well.

Author: Colton Bryant
*/

#include "stokesGrid.hpp"

//one-liners
/**
* Sets the wall velocity at y=yMax. Should be called during problem setup. \n
* Inputs: double val - the velocity you wish to set \n
* Result: uHigh set to val  
*/
void stokesGrid::setUppWallVelocity(const double val) {uHigh=val;} //upper wall velocity
/**
* Sets the wall velocity at y=yMin. Should be called during problem setup. \n 
* Inputs: double val - the velocity you wish to set \n
* Result: uLow set to val
*/
void stokesGrid::setLowWallVelocity(const double val) {uLow =val;} //lower wall velocity
/**
* Sets the relaxation parameter used in solving the momentum equation. Should be set during problem setup. \n
* Inputs: double val - the parameter to be used \n
* Result: omegaU set to val
*/
void stokesGrid::setOmegaU(const double val) {omegaU=val;} //momentum equation relaxation parameter
/**
* sets the relaxation parameter used in the pressure update. should be set during problem setup. \n
* inputs: double val - the parameter to be used \n
* result: omegaP set to val
*/
void stokesGrid::setOmegaP(const double val) {omegaP=val;} //pressure equation relaxation parameter
/**
* sets the tolerance for the SOR solver. should be set during problem setup. \n
* inputs: double val - the parameter to be used \n
* result: tol set to val
*/
void stokesGrid::setTol(const double val) {tol=val;} //tolerance for the iterative method

//grid point locations (xC = x-center location, xF = x-face location)
/**
* x coordinate of cell center. \n
* inputs: int i - index of cell \n
* return: x-coordinate of cell center
*/
double stokesGrid::xC(const int i) const {return mXMin+(i+0.5)*dx;}
/**
* y coordinate of cell center. \n
* inputs: int j - index of cell \n
* return: y-coordinate of cell center
*/
double stokesGrid::yC(const int j) const {return mYMin+(j+0.5)*dy;}
/**
* x coordinate of cell face. \n
* inputs: int i - index of cell \n
* return: x-coordinate of cell face
*/
double stokesGrid::xF(const int i) const {return mXMin+i*dx;}
/**
* y coordinate of cell face. \n
* inputs: int j - index of cell \n
* return: y-coordinate of cell face
*/
double stokesGrid::yF(const int j) const {return mXMin+j*dy;}


//functions for accessing flow data
/**
* x-velocity value. \n
* inputs: int i, int j - index of cell \n
* return: x-component of velocity field 
*/
double stokesGrid::u(const int i, const int j) const {return flowData[uI(i,j)];}
/**
* y-velocity value. \n
* inputs: int i, int j - index of cell \n
* return: y-component of velocity field 
*/
double stokesGrid::v(const int i, const int j) const {return flowData[vI(i,j)];}
/**
* pressure value. \n
* inputs: int i, int j - index of cell \n
* return: pressure in the cell
*/
double stokesGrid::p(const int i, const int j) const {return flowData[pI(i,j)];}

/**
* Converts i,j index to index in flowData array for x-velocity \n
* inputs: int i, int j - index of cell \n
* return: index of this value in flowData
*/
int stokesGrid::uI(const int i, const int j) const {return i*(mNy-1)+j;}
/**
* Converts i,j index to index in flowData array for y-velocity \n
* inputs: int i, int j - index of cell \n
* return: index of this value in flowData
*/
int stokesGrid::vI(const int i, const int j) const {return vStart+i*mNy+j;}
/**
* Converts i,j index to index in flowData array for pressure \n
* inputs: int i, int j - index of cell \n
* return: index of this value in flowData
*/
int stokesGrid::pI(const int i, const int j) const {return pStart+i*(mNy-1)+j;}


//functions to access body force data
/**
* Returns the x-component of the body force here\n
* inputs: int i, int j - index of cell \n
* return: x body force here
*/
double stokesGrid::fx(const int i, const int j) {return bodyXData[i*(mNy-1)+j];}
/**
* Returns the y-component of the body force here\n
* inputs: int i, int j - index of cell \n
* return: y body force here
*/
double stokesGrid::fy(const int i, const int j) {return bodyYData[i*mNy+j];}

//functions to set body force data 
/**
* Sets x body force in cell i,j \n
* inputs: int i, int j - index of cell; double val - value to set \n
* result: x-component of body force at index i,j set to val
*/
void stokesGrid::setFX(const int i, const int j, const double val) {bodyXData[i*(mNy-1)+j] = val;}
/**
* Sets y body force in cell i,j \n
* inputs: int i, int j - index of cell; double val - value to set \n
* result: y-component of body force at index i,j set to val
*/
void stokesGrid::setFY(const int i, const int j, const double val) {bodyYData[i*mNy+j] = val;}

//function to turn on verbose output
/**
* Turns on residual information for the Stokes solver \n
* inputs: none \n
* result: verbose set to true; any call of the "solve" function will write residual information to the terminal throughout the update
*/
void stokesGrid::turnOnVerbose() {verbose=true;}
//constructor
/**
*	Constructor for the Stokes grid. Specifies domain/grid size, allocates storage, and specifies default values for parameters. \n
*	Inputs: nx,ny - discrete grid size \n
*	xMin, xMax, yMin, yMax - bounds of the domain \n
*/
stokesGrid::stokesGrid(const int nx,
								const int ny,
								const double xMin,
								const double xMax,
								const double yMin,
								const double yMax)
{
	//set domain parameters
	mNx = nx;
	mNy = ny;
	mXMin = xMin;
	mYMin = yMin;
	mXMax = xMax;
	mYMax = yMax;

	//calculate dx,dy
	dx = (xMax-xMin)/(mNx-1);
	dy = (yMax-yMin)/(mNy-1);

	//set default values for tolerance and relaxation parameters. Each can be adjusted in main function.
	tol = 1.e-6;
	omegaU = 0.4;
	omegaP = 0.4;

	//default to no-slip/no-penetration conditions with wall velocities = 0
	uHigh = 0.;
	uLow  = 0.;

	//create data arrays
	flowData 	= new double[mNx*(mNy-1)+mNy*(mNx-1)+(mNx-1)*(mNy-1)]();
	bodyXData 	= new double[mNx*(mNy-1)]();
	bodyYData	= new double[mNy*(mNx-1)]();

	//set the starting locations for each variable in the flow data array
	uStart = 0;
	vStart = mNx*(mNy-1);
	pStart = vStart+mNy*(mNx-1);
	
	//default
	verbose = false;
}

//solver
/**
* Uses a modified SOR method to solve the stokes equations. \n
* Further details of the solve can be found in the mathematical documentation. 
*/
void stokesGrid::solve() {
	int count = 0;
	int i, j;
	int iP, iM;
	double here, west, east, nort, sout; //storage for data we pull down
	double tempRes; //temporary residual
	//prefactors
	double dxDy = dx/dy;
	double dyDx = dy/dx;
	double pMean;
	makeResidualsHuge();
	while(!isConverged() && count < MAXITS) {

		if(count%1000 == 0 && verbose) {
			cout << "On iteration " << count << " uRes = " << maxURes << ", vRes = " << maxVRes << ", pRes = " << maxPRes << endl;
		}
		zeroOutResiduals();
		//update x-velocity grid points
		for(i=0; i<mNx; ++i) {
			for(j=0; j<mNy-1; ++j) {
				if(i==mNx-1) {
					flowData[uI(i,j)] = u(0,j); //periodicity
				} else{
					iM = (i-1+mNx-1)%(mNx-1);
					iP = (i+1)%(mNx-1);
					here = u(i,j);	
					west = u(iM, j);
					east = u(iP, j);
					if(j==0) { //y=ymin
						nort = u(i,j+1);
						sout = 2.*uLow-u(i,j);
					} else if(j==mNy-2) { //y=ymax
						nort = 2.*uHigh-u(i,j);
						sout = u(i,j-1);
					} else { //other grid points
						nort = u(i,j+1);
						sout = u(i,j-1);
					}
					tempRes = dyDx * (west - 2.*here + east) + 
						dxDy * (sout - 2.*here + nort) - 
						dy * (p(i,j) - p(iM,j)) - 
						dx * dy * fx(i,j);

					flowData[uI(i,j)] += omegaU*tempRes;
					if( fabs(tempRes) > maxURes ) maxURes = fabs(tempRes);
				}
			}
		}
		//update y-velocity grid points
		for(i=0; i<mNx-1; ++i) {
			for(j=1; j<mNy-1; ++j) { //only update interior points b/c of no penetration condition
				iM = (i-1+mNx-1)%(mNx-1);
				iP = (i+1)%(mNx-1);
				here = v(i,j);
				west = v(iM,j);
				east = v(iP,j);
				nort = v(i,j+1);
				sout = v(i,j-1);
				tempRes = dyDx * (west - 2.*here + east) + 
					dxDy * (sout - 2.*here + nort) - 
					dx * (p(i,j) - p(i,j-1)) -
					dx * dy * fy(i,j);
				flowData[vI(i,j)] += omegaU*tempRes;
				if( fabs(tempRes) > maxVRes ) maxVRes = fabs(tempRes); 
			}
		}

		//update pressure
		pMean = 0.;
		for(i=0; i<mNx-1; ++i) {
			for(j=0; j<mNy-1; ++j) {
				tempRes = -(u(i+1,j)-u(i,j)) - dxDy*(v(i,j+1)-v(i,j));
				flowData[pI(i,j)] += omegaP*tempRes;
				pMean += flowData[pI(i,j)];
				if( fabs(tempRes) > maxPRes) maxPRes = fabs(tempRes);
			}
		}

		pMean /= (double) (mNx-1)*(mNy-1);
		for(i=0; i<mNx-1; ++i) {
			for(j=0; j<mNy-1; ++j) {
				flowData[pI(i,j)] -= pMean;
			}
		}

		count++;
	}
	if(count==MAXITS) {
		cout << "WARNING: Stokes solver failed to converge in " << MAXITS << " iterations.\n";
		cout << "Residuals: u = " << maxURes << ", v = " << maxVRes << ", p = " << maxPRes << endl;
		cout << "Exiting...\n";
		exit(EXIT_FAILURE);
	} else if(verbose) {
		cout << "Stokes solver converged in " << count << " iterations\n";
	}
}

//for forcing the solver into the loop
/**
* Sets residuals to large values to force us into the while loop in solve.
*/
void stokesGrid::makeResidualsHuge() {
	maxURes = 1e10;
	maxVRes = 1e10;
	maxPRes = 1e10;
}

//sets residuals to zero
/**
* Sets all maximum residuals back to zero. Done before each SOR update. 
*/
void stokesGrid::zeroOutResiduals() {
	maxURes = 0.;
	maxVRes = 0.;
	maxPRes = 0.;
}

//checks for convergence
/**
* Check convergence of the SOR solver. 
*/
bool stokesGrid::isConverged() {
	return maxURes < tol && maxVRes < tol && maxPRes < tol;
}

//interpolation routines
//generic bilinear interpolation routine
/**
* Generic bilinear interpolation routine. \n
* inputs: a11, a12, a21, a22 - data given on four corners of a unit square \n
* ifrac, jfrac - values between 0 and 1 specifying the location to interpolate to. \n
* e.g. ifrac = 0.5, jfrac = 0.75 will interpolate to the point (0.5, 0.75) \n
* return: interpolated value at the given location
*/
double stokesGrid::bilin(const double a11, const double a12, const double a21, const double a22, 
	const double ifrac, const double jfrac) const 
{
	return (1.-ifrac)*(1.-jfrac)*a11 + 
		(1.-ifrac)*jfrac*a12 + 
		ifrac*(1.-jfrac)*a21 + 
		ifrac*jfrac*a22;
}
//flow variable interpolation
/**
* Bilinear interpolation for the x-velocity. \n
* Converts a given point (x,y) to an x-velocity cell of the staggered grid. Then uses bilin to interpolate. Note: 
* this function accounts for the BCs. \n
* inputs: double x,y - location in the domain \n
* return: interpolated value of u(x,y) \n
*/
double stokesGrid::uB(const double x, const double y) const {
	int i, j, iP;
	double ifrac, jfrac;
	//find the grid cell containing this point
	i = (int) floor(((x-mXMin)/dx));
	j = (int) floor(((y-(mYMin+dy/2))/dy));
	
	ifrac = (x-xF(i))/dx;
	jfrac = (y-yC(j))/dy;

	//set iP enforcing periodicty
	iP = i < mNx-1 ? i+1 : 1;

	//cases for j
	if(j == -1) { //y \in (yMin, yMin+dy/2)
		return bilin(2.*uLow-u(i,0), u(i,j+1), 2.*uLow-u(iP,0), u(iP,j+1), ifrac, jfrac);
	} else if(j == mNy-2) { //y \in(yMax-dy/2, yMax)
		return bilin(u(i,j), 2.*uHigh-u(i,j), u(iP,j), 2.*uHigh-u(iP,j), ifrac, jfrac);
	} else {
		return bilin(u(i,j), u(i,j+1), u(i+1,j), u(i+1,j+1), ifrac, jfrac);
	}
}
/**
* Bilinear interpolation for the y-velocity. \n
* Converts a given point (x,y) to an x-velocity cell of the staggered grid. Then uses bilin to interpolate. Note: 
* this function accounts for the BCs. \n
* inputs: double x,y - location in the domain \n
* return: interpolated value of v(x,y) \n
*/
double stokesGrid::vB(const double x, const double y) const {
	int i, j, iP;
	double ifrac, jfrac;
	i = (int) floor(((x-(mXMin+dx/2))/dx));
	j = (int) floor(((y-mYMin)/dy));
	
	ifrac = (x-xC(i))/dx;
	jfrac = (y-yF(j))/dy;

	iP = i < mNx-2 ? i+1 : 1;
	i = i < 0 ? mNx-2 : i;

	return bilin(v(i,j), v(i,j+1), v(iP,j), v(iP, j+1), ifrac, jfrac);
	
}
/**
* Bilinear interpolation for the pressure. \n
* Converts a given point (x,y) to an x-velocity cell of the staggered grid. Then uses bilin to interpolate. Note:
* due to lack of pressure data at y=ymin, y=ymax, pressure cannot be interpolated to the upper/lower walls. \n
* inputs: double x,y - location in the domain \n
* return: interpolated value of p(x,y) \n
*/
double stokesGrid::pB(const double x, const double y) const {
	int i, j, iP;
	double ifrac, jfrac;
	i = (int) floor(((x-(mXMin+dx/2))/dx));
	j = (int) floor(((y-(mYMin+dy/2))/dy));
	ifrac = (x-xC(i))/dx;
	jfrac = (y-yC(j))/dy;
	
	iP = i < mNx-2 ? i+1 : 1;
	i = i < 0 ? mNx-2 : i;
	
	//note pressure cannot be interpolated up to the walls. This is due to lack of BCs. If you need pressure there you'll
	//need to come up with something
	return bilin(p(i,j), p(i,j+1), p(iP,j), p(iP,j+1), ifrac, jfrac);

}



//data dumps
/**
* Writes velocity and pressure data to binary files. \n
* Creates 3 output files, one for x-velocity, one for y-velocity, one for pressure \n
* Each file contains the following data: \n
* integer xLen: the number of grid points in x for this array \n
* integer yLen: the number of grid points in y for this array \n
* double data: an xLen*yLen length array of the data \n
* double coords: an xLen+yLen length array containing the x then y coordinates of nodes where this data was stored in the staggered grid. \n
* \n
* inputs: string uOut, vOut, pOut - filenames for the data to be written to.
* result: data written as described to these files in binary form
*/
void stokesGrid::dumpFlowData(string uOut, string vOut, string pOut) {
	FILE* uFile;
	FILE* vFile;
	FILE* pFile;
	uFile = fopen(uOut.c_str(), "w");
	vFile = fopen(vOut.c_str(), "w");
	pFile = fopen(pOut.c_str(), "w");
	
	int uxLen, uyLen, vxLen, vyLen, pxLen, pyLen;
	uxLen = mNx;
	uyLen = mNy-1;
	vxLen = mNx-1;
	vyLen = mNy;
	pxLen = mNx-1;
	pyLen = mNy-1;
	
	//write array sizes
	fwrite(&uxLen, sizeof(int), 1, uFile);
	fwrite(&uyLen, sizeof(int), 1, uFile);
	fwrite(&vxLen, sizeof(int), 1, vFile);
	fwrite(&vyLen, sizeof(int), 1, vFile);
	fwrite(&pxLen, sizeof(int), 1, pFile);
	fwrite(&pyLen, sizeof(int), 1, pFile);

	//write data
	fwrite(&(flowData[uStart]), sizeof(double), uxLen*uyLen, uFile);
	fwrite(&(flowData[vStart]), sizeof(double), vxLen*vyLen, vFile);
	fwrite(&(flowData[pStart]), sizeof(double), pxLen*pyLen, pFile);

	//write coordinates
	double* uCoords = new double[uxLen+uyLen];
	double* vCoords = new double[vxLen+vyLen];
	double* pCoords = new double[pxLen+pyLen];
	
	int i,j;
	for(i=0; i<uxLen; ++i) uCoords[i] = xF(i);
	for(j=0; j<uyLen; ++j) uCoords[uxLen+j] = yC(j);
	for(i=0; i<vxLen; ++i) vCoords[i] = xC(i);
	for(j=0; j<vyLen; ++j) vCoords[vxLen+j] = yF(j);
	for(i=0; i<pxLen; ++i) pCoords[i] = xC(i);
	for(j=0; j<pyLen; ++j) pCoords[pxLen+j] = yC(j);
	
	fwrite(uCoords, sizeof(double), uxLen+uyLen, uFile);
	fwrite(vCoords, sizeof(double), vxLen+vyLen, vFile);
	fwrite(pCoords, sizeof(double), pxLen+pyLen, pFile);

	delete[] uCoords;
	delete[] vCoords;
	delete[] pCoords;

	fclose(uFile);
	fclose(vFile);
	fclose(pFile);
}
/**
* Writes the x-coordinate of the body force. See dumpFlowData for format information.
*/
void stokesGrid::dumpFX(string fxOut) {
	FILE* fxFile;
	fxFile = fopen(fxOut.c_str(), "w");
	
	int xLen = mNx;
	int yLen = mNy-1;

	//write size
	fwrite(&xLen, sizeof(int), 1, fxFile);
	fwrite(&yLen, sizeof(int), 1, fxFile);
	
	//write data
	fwrite(bodyXData, sizeof(double), xLen*yLen, fxFile);

	double* coords = new double[xLen+yLen];
	for(int i=0; i<xLen; ++i) coords[i] = xF(i);
	for(int j=0; j<yLen; ++j) coords[xLen+j] = yC(j);
	
	fwrite(coords, sizeof(double), xLen+yLen, fxFile);
	delete[] coords;
	fclose(fxFile);
}

/**
* Writes the y-coordinate of the body force. See dumpFlowData for format information.
*/
void stokesGrid::dumpFY(string fyOut) {
	FILE* fyFile;
	fyFile = fopen(fyOut.c_str(), "w");
	
	int xLen = mNx-1;
	int yLen = mNy;

	//write size
	fwrite(&xLen, sizeof(int), 1, fyFile);
	fwrite(&yLen, sizeof(int), 1, fyFile);
	
	//write data
	fwrite(bodyYData, sizeof(double), xLen*yLen, fyFile);

	double* coords = new double[xLen+yLen];
	for(int i=0; i<xLen; ++i) coords[i] = xC(i);
	for(int j=0; j<yLen; ++j) coords[xLen+j] = yF(j);
	
	fwrite(coords, sizeof(double), xLen+yLen, fyFile);
	delete[] coords;
	fclose(fyFile);
}


//destructor
/**
* Deallocates storage. 
*/
stokesGrid::~stokesGrid(void) {
	if(flowData) delete[] flowData;
	if(bodyXData) delete[] bodyXData;
	if(bodyYData) delete[] bodyYData;
}









