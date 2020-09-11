/*
Definitions for the stokes grid class.

Notes: It is assumed that boundary conditions are periodic in x and uses no slip/no penetration conditions in y. Non-zero wall velocites can be
specified at y=ymin and y=ymax as well.

Author: Colton Bryant
*/

#include "stokesGrid.hpp"

//one-liners
void stokesGrid::setUppWallVelocity(const double val) {uHigh=val;} //upper wall velocity
void stokesGrid::setLowWallVelocity(const double val) {uLow =val;} //lower wall velocity

void stokesGrid::setOmegaU(const double val) {omegaU=val;} //momentum equation relaxation parameter
void stokesGrid::setOmegaP(const double val) {omegaP=val;} //pressure equation relaxation parameter
void stokesGrid::setTol(const double val) {tol=val;} //tolerance for the iterative method

//grid point locations (xC = x-center location, xF = x-face location)
double stokesGrid::xC(const int i) const {return mXMin+(i+0.5)*dx;}
double stokesGrid::yC(const int j) const {return mYMin+(j+0.5)*dy;}
double stokesGrid::xF(const int i) const {return mXMin+i*dx;}
double stokesGrid::yF(const int j) const {return mXMin+j*dy;}


//functions for accessing flow data
double stokesGrid::u(const int i, const int j) const {return flowData[uI(i,j)];}
double stokesGrid::v(const int i, const int j) const {return flowData[vI(i,j)];}
double stokesGrid::p(const int i, const int j) const {return flowData[pI(i,j)];}

int stokesGrid::uI(const int i, const int j) const {return i*(mNy-1)+j;}
int stokesGrid::vI(const int i, const int j) const {return vStart+i*mNy+j;}
int stokesGrid::pI(const int i, const int j) const {return pStart+i*(mNy-1)+j;}


//functions to access body force data
double stokesGrid::fx(const int i, const int j) {return bodyXData[i*(mNy-1)+j];}
double stokesGrid::fy(const int i, const int j) {return bodyYData[i*mNy+j];}

//functions to set body force data 
void stokesGrid::setFX(const int i, const int j, const double val) {bodyXData[i*(mNy-1)+j] = val;}
void stokesGrid::setFY(const int i, const int j, const double val) {bodyYData[i*mNy+j] = val;}

//function to turn on verbose output
void stokesGrid::turnOnVerbose() {verbose=true;}
//constructor
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
/*
Uses a modified SOR method to solve the stokes equations.
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
void stokesGrid::makeResidualsHuge() {
	maxURes = 1e10;
	maxVRes = 1e10;
	maxPRes = 1e10;
}

//sets residuals to zero
void stokesGrid::zeroOutResiduals() {
	maxURes = 0.;
	maxVRes = 0.;
	maxPRes = 0.;
}

//checks for convergence
bool stokesGrid::isConverged() {
	return maxURes < tol && maxVRes < tol && maxPRes < tol;
}

//interpolation routines
//generic bilinear interpolation routine
double stokesGrid::bilin(const double a11, const double a12, const double a21, const double a22, 
	const double ifrac, const double jfrac) const 
{
	return (1.-ifrac)*(1.-jfrac)*a11 + 
		(1.-ifrac)*jfrac*a12 + 
		ifrac*(1.-jfrac)*a21 + 
		ifrac*jfrac*a22;
}
//flow variable interpolation
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
stokesGrid::~stokesGrid(void) {
	if(flowData) delete[] flowData;
	if(bodyXData) delete[] bodyXData;
	if(bodyYData) delete[] bodyYData;
}









