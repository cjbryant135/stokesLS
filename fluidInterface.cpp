/*
Definitions for the fluid interface class.

Author: Colton Bryant
*/

#include "fluidInterface.hpp"

//one-liners
/**
* Set the surface tension parameter. \n
* inputs: double sigma - surface tension parameter \n
* result: mSigma set to sigma. 
*/
void fluidInterface::setST(const double sigma) {mSigma = sigma;} //set surface tension
/**
* Set the width of the discrete delta function used in computing surface tension as a body force. \n
* inputs: double eps - width to use \n
* result: mEps set to eps. 
*/
void fluidInterface::setDeltaWidth(const double eps) {mEps = eps;}; //width of the discrete delta used in surface tension calculations

//getters for domain size
/**
* Get domain bounds. Used this for debugging sometimes \n
* inputs: none \n
* return: xMin 
*/
double fluidInterface::getXMin() {return grid->X(0);} 
/**
* Get domain bounds. Used this for debugging sometimes \n
* inputs: none \n
* return: yMin 
*/
double fluidInterface::getYMin() {return grid->Y(0);}
/**
* Get domain bounds. Used this for debugging sometimes \n
* inputs: none \n
* return: xMax 
*/
double fluidInterface::getXMax() {return grid->X(mNx-1);}
/**
* Get domain bounds. Used this for debugging sometimes \n
* inputs: none \n
* return: yMax 
*/
double fluidInterface::getYMax() {return grid->Y(mNy-1);}

//get the index for a given slice in the level set code
/**
* Returns the "slice index" corresponding to a valid named slice. Slices are named by strings in the constructor. \n
* input: string slice - name of a slice of the level set data \n
* return: k-index corresponding to that slice in the level set code. Note: An invalid string will throw an error.
*/
int fluidInterface::nameToIndx(const string slice) {return dIndx.at(slice);} //use "at" so an error is thrown for invalid slice names
//dump a named slice of the level set code. See constructor for available names
/**
* Writes a named slice of the level set data to file. \n
* input: string slice - name of a slice of the level set data, string filename - filename where data is to be written \n
* result: level set data for the given slice is written to file. See LevelSet/src/um2io.cpp for details of the output file.
*/
void fluidInterface::dumpSlice(const string slice, const string filename) {grid->WriteBinary(filename, nameToIndx(slice));}

//constructor
/**
* Creates the level set grid for our model of the fluid/fluid interface. A map is used to name each "slice" of the level set data. \n
* Inputs: nx, ny - discrete grid size \n
* xMin, xMax, yMin, yMax - domain bounds \n
* sigma - surface tension \n
*/
fluidInterface::fluidInterface(	const int nx,
											const int ny,
											const double xMin,
											const double xMax,
											const double yMin,
											const double yMax,
											const double sigma)
{
	mNx = nx;
	mNy = ny;
	//create the double slice index map: The level set code allows derivatives and 
	//other data to be computed over the entire mesh and stored in a slice. Here we
	//map these slices to strings for easy access later.
	dIndx["phi"] 		= 0; //level set function phi(x,y)
	dIndx["phi_x"] 	= 1; //phi_x (centered differencing)
	dIndx["phi_y"] 	= 2; //phi_y (centered differencing)
	dIndx["phi_xx"] 	= 3; //phi_xx (centered differencing)
	dIndx["phi_xy"] 	= 4; //phi_xy (centered differencing)
	dIndx["phi_yy"] 	= 5; //phi_yy (centered differencing)
	dIndx["kappa"] 	= 6; //curvature 
	dIndx["surf_x"] 	= 7; //x component of surface tension force
	dIndx["surf_y"] 	= 8; //y component of surface tension force
	dIndx["speed"]    = 9; //extended interface speed
	dIndx["phi_cpy"]	= 10; //copy of phi for the velocity extension
	dIndx["phi_xp"] 	= 11; //phi_x (forward differencing)
	dIndx["phi_yp"] 	= 12; //phi_y (forward differencing)
	dIndx["phi_xm"] 	= 13; //phi_x (backward differencing)
	dIndx["phi_ym"] 	= 14; //phi_y (backward differencing)
	dIndx["phi_ng"]	= 15; // \| \grad \phi \| 
	//create the integer slice index map: Similar to above but for integer data
	iIndx["mask"] = 0; //locations of initialized speed data for velocity extension
	iIndx["ikus"] = 1; //update status array for the fast marching method
	iIndx["ikhi"] = 2; //integer storage for fast marching
	
	//set total number of slices
	dSlices = 16; //double slices
	iSlices = 3; //integer slices
	
	//set surface tension
	mSigma = sigma;

	//create the grid
	grid = new UniformMesh2D(Bounds, nx, ny, dSlices, iSlices, xMin, xMax, yMin, yMax, *(new UM2_LinearBdry(2,2,2,2)));

	//set default eps to be 1.5 dx
	mEps = 1.5*grid->dx;
}

//initialize the level set function
/**
* Sets initial level set function data using an InitialFunc from the level set code \n
* Inputs: InitialFunc f - object specifying interface shape from the level set code. Other shapes found in 
*								LevelSet/src/initfuncs/ \n
* Result: the "phi" slice of level set data is set to a signed distance function with zero level set on the given shape. 
*/
void fluidInterface::initialize(const InitialFunc& f) {
	grid->SetValues(f, dIndx.at("phi")); //call the set values function to initialize phi as a signed distance function
}

/**
* Computes the surface tension body force using a discrete delta function of width 2*eps. Note: This is computed on the 
* nodes of the level set grid and should later be interpolated to the staggered grid for use with the Stokes solver.
*/
void fluidInterface::computeSurfaceTension() {
	int i,j;
	int sX = dIndx.at("surf_x");
	int sY = dIndx.at("surf_y");
	//compute derivatives of phi 
	grid->Dx_zero(dIndx.at("phi"), dIndx.at("phi_x"));
	grid->Dy_zero(dIndx.at("phi"), dIndx.at("phi_y"));
	grid->Dxx_zero(dIndx.at("phi"), dIndx.at("phi_xx"));
	grid->Dxy_zero(dIndx.at("phi"), dIndx.at("phi_xy"));
	grid->Dyy_zero(dIndx.at("phi"), dIndx.at("phi_yy"));
	//compute curvature
	grid->Curvature(dIndx.at("phi_x"), dIndx.at("phi_y"), 
		dIndx.at("phi_xx"), dIndx.at("phi_xy"), dIndx.at("phi_yy"), dIndx.at("kappa"));
	//compute surface tension
	for(i=0; i<grid->maxi; ++i) {
		for(j=0; j<grid->maxj; ++j) {
			//x component
			grid->data_(i,j,sX) = fabs(grid->data_(i,j,dIndx.at("phi"))) < mEps ? 
				mSigma*grid->data_(i,j,dIndx.at("phi_x"))*grid->data_(i,j,dIndx.at("kappa"))*(1.+cos(M_PI*grid->data_(i,j,dIndx.at("phi"))/mEps))/2./mEps : 0.;
			//y component
			grid->data_(i,j,sY) = fabs(grid->data_(i,j,dIndx.at("phi"))) < mEps ? 
				mSigma*grid->data_(i,j,dIndx.at("phi_y"))*grid->data_(i,j,dIndx.at("kappa"))*(1.+cos(M_PI*grid->data_(i,j,dIndx.at("phi"))/mEps))/2./mEps : 0.;
		}
	}
	//apply BCs
	grid->bc->Apply(sX);
	grid->bc->Apply(sY);
}

/** 
* Computes surface tension and interpolates it onto the staggered grid for use by the stokes solver \n
* Inputs: sGrid - a stokesGrid object \n
*/
void fluidInterface::computeSurfaceTension(stokesGrid* sGrid) {
	//compute surface tension on the level set grid
	computeSurfaceTension();
	//populate arrays on the staggered grid
	int i,j;
	for(i=0; i<mNx; ++i) {
		for(j=0; j<mNy; ++j) {
			if( j < mNy-1 ) {
				//set x component of body force
				sGrid->setFX(i,j,(grid->data_(i,j,dIndx.at("surf_x"))+grid->data_(i,j+1,dIndx.at("surf_x")))/2.);
			}
			if( i < mNx-1 ) {
				//set y component of body force
				sGrid->setFY(i,j,(grid->data_(i,j,dIndx.at("surf_y"))+grid->data_(i+1,j,dIndx.at("surf_y")))/2.);
			}
		}
	}
}

/**
* Velocity extension routine. Extrapolates normal velocity at the interface to the rest of the domain by following grad phi.
* Note: This function assumes x and y derivatives of phi have already been computed and stored in "phi_x" and "phi_y"
* slices. This happens automatically if surface tension has been computed earlier in the timestep. \n
* Inputs: stokesGrid sGrid - a stokesGrid object whose velocity data will be used to compute the interface speed\n
* Result: the "speed" slice of the level set data contains the extensional velocity F and "phi" can now be updated.
*/
void fluidInterface::setSpeed(const stokesGrid* sGrid) {
	//compute norm grad phi
	grid->NormGrad(dIndx.at("phi_x"), dIndx.at("phi_y"), dIndx.at("phi_ng"));

	//some integers
	int i, j, ii, jj;
	int cross;

	//zero out mask
	for(i=0; i<grid->maxi; ++i) {
		for(j=0; j<grid->maxj; ++j) {
			grid->idata_(i,j,iIndx.at("mask")) = 0;
		}
	}
	
	//Initialize the speed near the interface

	//initialze storage for a bicubic
	levelset::Bicubic p;
	//some storage
	double ax, ay;
	double dist;
	bool clean;
	double uTemp, vTemp;
	//loop over voxels of the domain
	for(i=0; i<grid->maxi-1; ++i) {
		for(j=0; j<grid->maxj-1; ++j) {
			//check for a sign change in phi ==> the interface crosses this grid cell
			cross = 0;
			cross += grid->data_(i,j,dIndx.at("phi")) > 0. ? 1 : 0;
			cross += grid->data_(i+1,j,dIndx.at("phi")) > 0. ? 2 : 0;
			cross += grid->data_(i,j+1,dIndx.at("phi")) > 0. ? 4 : 0;
			cross += grid->data_(i+1,j+1,dIndx.at("phi")) > 0. ? 8 : 0;
			if(cross > 0 && cross < 15) { //sign change occured
				//build a bicubic fit to phi over this cell
				p.BuildwDeriv( 
					grid->data_(i-1,j-1,dIndx.at("phi")), grid->data_(i,j-1,dIndx.at("phi")), grid->data_(i+1,j-1,dIndx.at("phi")), grid->data_(i+2,j-1,dIndx.at("phi")), 
					grid->data_(i-1,j,dIndx.at("phi")), grid->data_(i,j,dIndx.at("phi")), grid->data_(i+1,j,dIndx.at("phi")), grid->data_(i+2,j,dIndx.at("phi")), 
					grid->data_(i-1,j+1,dIndx.at("phi")), grid->data_(i,j+1,dIndx.at("phi")), grid->data_(i+1,j+1,dIndx.at("phi")), grid->data_(i+2,j+1,dIndx.at("phi")), 
					grid->data_(i-1,j+2,dIndx.at("phi")), grid->data_(i,j+2,dIndx.at("phi")), grid->data_(i+1,j+2,dIndx.at("phi")), grid->data_(i+2,j+2,dIndx.at("phi")),
					grid->dx, grid->dy); 
				//loop over corners
				for(ii=i; ii <= i+1; ++ii) {
					for(jj=j; jj <= j+1; ++jj) {
						//compute distance to the interface using the modified newton
						dist = distToInt(&p, grid->X(ii), grid->Y(jj), grid->X(i), grid->X(i+1), 
							grid->Y(j), grid->Y(j+1), ax, ay, clean);
						//check if this converged 
						if(clean) {
							//set mask to one here
							grid->idata_(ii,jj,iIndx.at("mask")) = 1;
							//set speed 
							grid->data_(ii,jj,dIndx.at("speed")) = (
								sGrid->uB(ax, ay)*grid->Interp(ax, ay, dIndx.at("phi_x")) + 
								sGrid->vB(ax, ay)*grid->Interp(ax, ay, dIndx.at("phi_y")) )/
								grid->Interp(ax, ay, dIndx.at("phi_ng") );
						}
					}
				}
				
			}
		} 
	} //end of voxel loop

	//make a copy of phi (note this copy may be destroyed after velocity extension)
	grid->CopyWorkGrid(dIndx.at("phi_cpy"), dIndx.at("phi"));
	
	//call the velocity extension routine (NOTE: There are several versions of this routine in the code. 
	//Here the version used by Chimera code is called even though it may not be necessary. Other versions 
	//can be found in LevelSet/src/um2advect.cpp
	grid->ExtendVelocity(dIndx.at("phi"), dIndx.at("phi_cpy"), dIndx.at("speed"), 
		iIndx.at("ikus"), iIndx.at("ikhi"), iIndx.at("mask"), false, false);
	grid->bc->Apply(dIndx.at("speed"));

}
/**
* Uses Newtons method to find the closest point on the interface to a given grid point and the distance to it \n
* Inputs: Bicubic* fit - pointer to a bicubic fit of the level set data near the interface \n
	xi, yj - Location of the grid point we are computing the distance to \n
	xMin, xMax, yMin, yMax - Corners of the grid cell of interest \n
	xC, yC - Location of the root found by the Newton's method (passed by reference for use later) \n
	clean - set to true if the newton's method converges in maxIts iterations AND the root found is inside the current grid cell \n
* Return: distance to the interface
*/
double fluidInterface::distToInt(levelset::Bicubic* fit, const double xi, const double yj, const double xMin, const double xMax, 
	const double yMin, const double yMax, double &xC, double &yC, bool &clean) 
{

	//initialize some storage 
	double p, pX, pY, pXX, pYY, pXY, d1, d2, xDiff, yDiff;
	clean = false;

	int maxIts = 100; //max number of iterations for the Newton Solver
	double tol = 1e-5; //tolerance for the Newton solver

	//set an initial guess
	xC = xi;
	yC = yj;


	int count = 0;

	double deltaX, deltaY;
	deltaX = 1000.;
	deltaY = 1000.;

	while( (deltaX > tol || deltaY > tol) && count < maxIts ) {
		//compute terms needed for the Jacobian
		p 		= fit->F(xC-xMin,yC-yMin);
		pX 	= fit->Dx(xC-xMin, yC-yMin);
		pY    = fit->Dy(xC-xMin, yC-yMin);
		pXX   = fit->Dxx(xC-xMin, yC-yMin);
		pYY   = fit->Dyy(xC-xMin, yC-yMin);
		pXY   = fit->Dxy(xC-xMin, yC-yMin);
	
		//additional variables
		xDiff = xC-xi;
		yDiff = yC-yj;

		d1 = pYY*xDiff - pXY*yDiff - pX;
		d2 = pY + pXY*xDiff - pXX*yDiff;

		//compute the Newton update (we are direclty using the 2x2 inverse formula here)
		deltaX = (p*d1 - pY*pY*xDiff + pX*pY*yDiff)/(pX*d1-pY*d2);
		deltaY = (-p*d2 + pX*pY*xDiff - pX*pX*yDiff)/(pX*d1-pY*d2);

		//update 
		xC -= deltaX;
		yC -= deltaY;
		
		//increase iterations
		count++;
	}
	
	clean = (xC >= xMin && xC <= xMax && yC >= yMin && yC <= yMax && count < maxIts);

	
	return sqrt((xC-xi)*(xC-xi) + (yC-yj)*(yC-yj));

}

/**
* Use upwinding to advance the interface according to speed stored in the "speed" slice\n
* Inputs: dt - the timestep to use \n
* result: "phi" data updated according to the upwind scheme
*/
void fluidInterface::upwindAdvance(const double dt) {
	int i,j;
	double upwindGrad;
	double tempx, tempy, F;
	//compute forward/backward differences
	grid->Dx_minus(dIndx.at("phi"), dIndx.at("phi_xm"));
	grid->Dx_plus(dIndx.at("phi"), dIndx.at("phi_xp"));
	grid->Dy_minus(dIndx.at("phi"), dIndx.at("phi_ym"));
	grid->Dy_plus(dIndx.at("phi"), dIndx.at("phi_yp"));
	
	for(i=0; i<grid->maxi; ++i) {
		for(j=0; j<grid->maxj; ++j) {
			F 		= grid->data_(i,j,dIndx.at("speed"));
			tempx = std::max(sgn(F)*grid->data_(i,j,dIndx.at("phi_xm")),-sgn(F)*grid->data_(i,j,dIndx.at("phi_xp")));
			tempy = std::max(sgn(F)*grid->data_(i,j,dIndx.at("phi_ym")),-sgn(F)*grid->data_(i,j,dIndx.at("phi_yp")));
			
			tempx = tempx > 0. ? tempx : 0.;
			tempy = tempy > 0. ? tempy : 0.;

			upwindGrad = F*sqrt(tempx*tempx + tempy*tempy);

			grid->data_(i,j,dIndx.at("phi")) -= dt*upwindGrad;
		}
	}
	
	grid->bc->Apply(dIndx.at("phi"));

}

/**
* Deallocates the level set data
*/
fluidInterface::~fluidInterface() {
	if(grid) delete grid;
}

