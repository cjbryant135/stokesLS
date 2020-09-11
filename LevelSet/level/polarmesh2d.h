//
//  polarmesh2d.h
//  LevelSet
//
//  Created by David Chopp on 7/23/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#ifndef polarmesh2d_h
#define polarmesh2d_h

#ifndef NO_GRAPHICS
#include "plotwindow2d.h"
#include "plotwindow3d.h"
#endif
#include "initialfunc.h"
#include "genbicubic.h"
#include "datavec.h"
#include "utility.h"
#include "interface.h"
#include "numtrait.h"
#include <string>
#include <math.h>
#include <limits.h>
//#include "fluidbc.h"
#ifdef LEVEL_DEBUG
#include "debug.h"
#endif

namespace levelset {
    
#ifndef MESH_ENUMS_DEFINED
#define MESH_ENUMS_DEFINED
    enum STEPSIZE {StepSize};
    enum BOUNDS {Bounds};
    enum Direction {DLeft = 1, DDown = 2, DRight = 4, DUp = 8};
    enum {Upwind = -1, Downwind = 1};
#endif
    
    class PM2_HeapElement;
    class PolarMesh2D_Heap;
    class PM2_Boundary;
    
    class PolarMesh2D {
    public:
        
        double      *thedata;            // gridpoint array
        int         datastart;
        int         maxi, maxj, maxk;    // dimensions of data
        int         tmaxi, tmaxj;        // total dimensions of data
        double      *r, *dr;
        double      *tr, *tdr;
        double      dtheta;              // spatial step sizes
        int         *theidata;           // integer arrays for grid indexing heap
        int         nextw;
        double      *fmm_coeff[4][4];
        bool        *fmm_coi;
        double      *fmm_val[4][4];
        double      *fmm_ext[4][3];
        PM2_Boundary* bc;
        
        inline int Index(const int i, const int j, const int k=0) const
        {return(k*tmaxi*tmaxj+i*tmaxj+j);}
        inline int Index00(const int i, const int j, const int k=0) const
        {return(k*tmaxi*tmaxj+i*tmaxj+j+datastart);}
        
    public:
        
        PolarMesh2D(STEPSIZE, const int m, const int n,
                    const int sd, const int si,
                    const double deltar,
                    const double rmin, const double rstretch,
                    PM2_Boundary& thebc);
        
        PolarMesh2D(BOUNDS, const int m, const int n,
                    const int sd, const int si,
                    const double rmin, const double rmax,
                    const double rstretch,
                    PM2_Boundary& thebc);
        
        PolarMesh2D(void) {};
        
        virtual ~PolarMesh2D(void);
        
        // Initialization routines
        
        void SetStepSize(const int m, const int n,
                         const int sd, const int si,
                         const double deltar,
                         const double r0, PM2_Boundary& thebc);
        
        void SetBoundsSize(const int m, const int n,
                           const int sd, const int si,
                           const double r0, const double r1,
                           PM2_Boundary& thebc);
        
        void BuildWorkGrid(int num);
        void CopyWorkGrid(int to, int from);
		  void InitFMMCoeffs(void);
        
        void SetValuesXY(const InitialFunc& f, const int wi, const double mult = 1.);
        void SetValuesRT(const InitialFunc& f, const int wi, const double mult = 1.);
        void SetValuesXY(double (*f)(double,double), const int wi, const double mult = 1.);
        void SetValuesRT(double (*f)(double,double), const int wi, const double mult = 1.);
        void SetValues(const double* somedata, const int wi, const double mult = 1.);
        void SetValues(const double c, const int wi);

        // Data Access
#ifdef USE_INLINES
        inline double& data_(const int i, const int j, const int k)
        {int t=Index00(i,j,k); return(thedata[t]);}
        inline double data_(const int i, const int j, const int k) const
        {return(thedata[Index00(i,j,k)]);}
        inline int& idata_(const int i, const int j, const int k)
        {return(theidata[Index00(i,j,k)]);}
        inline int idata_(const int i, const int j, const int k) const
        {return(theidata[Index00(i,j,k)]);}
#define mdata_(M,I,J,K) data_(I,J,K)
#define midata_(M,I,J,K) idata_(I,J,K)
#else
#define data_(I,J,K) thedata[(K)*tmaxi*tmaxj+(I)*tmaxj+(J)+datastart]
#define mdata_(M,I,J,K) thedata[(K)*(M->tmaxi)*(M->tmaxj)+(I)*(M->tmaxj)+(J)+(M->datastart)]
#define idata_(I,J,K) theidata[(K)*tmaxi*tmaxj+(I)*tmaxj+(J)+datastart]
#define midata_(M,I,J,K) theidata[(K)*(M->tmaxi)*(M->tmaxj)+(I)*(M->tmaxj)+(J)+(M->datastart)]
#endif
        
        // Information retrieval
        
        double X(const int i, const int j) const;
        double Y(const int i, const int j) const;
        double X(const double i, const double j) const;
        double Y(const double i, const double j) const;
        int I(const double x, const double y) const;
        inline int J(const double x, const double y) const {double theta=atan2(y, x); 
		  	theta += theta <= 0. ? 2.*M_PI : 0. ; return int(floor(theta/dtheta));}
        double Id(const double x, const double y) const;
        inline double Jd(const double x, const double y) const {
		  		double theta=atan2(y, x); 
				//std::cout << "theta org = " << theta << std::endl;
				theta += theta <= 0. ? 2.*M_PI : 0. ;
		  		//std::cout << "theta shi = " << theta << std::endl;
				return theta/dtheta;
			}

        inline double Rho(const int i) const {return r[i];}
        inline double Theta(const int j) const {return j*dtheta;}
        inline double Rho(const double i) const {return r[int(i)]+(i-int(i))*dr[int(i)];}
        inline double Theta(const double j) const {return j*dtheta;}
        int I(const double rho) const;
        inline int J(const double theta) const {return int(floor(theta/dtheta));}
        double Id(const double rho) const;
        inline double Jd(const double theta) const {return theta/dtheta;}
        
        void LocToIndexXY(const double x, const double y, int& i, int& j,
                          double& ifrac, double& jfrac, const char round=0) const;
        void LocToIndexRT(const double rho, const double theta, int& i, int& j,
                          double& ifrac, double& jfrac, const char round=0) const;

        double InterpXY(const double& x, const double& y, const int k) const;
        double InterpRT(const double& rho, const double& theta, const int k) const;
        double Interp2XY(const double& x, const double& y, const int k) const;
        double Interp2RT(const double& rho, const double& theta, const int k) const;
        double Interp2RT(const int i, const int j, const int k, const double& ifrac,
                         const double& jfrac) const;
        GeneralBicubic Interp2(const int i, const int j, const int k) const;

        // I/O routines
        
        void ReadBinaryXY(const std::string fname, const int kstart = 0);
        void ReadBinaryXY(const char* fname, const int kstart = 0);
        void WriteBinaryXY(const std::string fname, const int kstart = -1,
                         const int kend = -1, const char* comment = NULL);
        void WriteBinaryXY(const char* fname, const int kstart = -1,
                         const int kend = -1, const char* comment = NULL);
        void WriteMatlabCoordinatesXY(const char* fname);

        void ReadBinaryRT(const std::string fname, const int kstart = 0);
        void ReadBinaryRT(const char* fname, const int kstart = 0);
        void WriteBinaryRT(const std::string fname, const int kstart = -1,
                           const int kend = -1, const char* comment = NULL);
        void WriteBinaryRT(const char* fname, const int kstart = -1,
                           const int kend = -1, const char* comment = NULL);
        void WriteMatlabCoordinatesRT(const char* fname);
        
        // Whole mesh methods
        
        void Dr_zero(const int k, const int kdr);
        void Dr_minus(const int k, const int kdr);
        void Dr_plus(const int k, const int kdr);
        void Dr_upwind(const int v, const int k, const int kdr);
        void Drr_zero(const int k, const int kdr);
        
        void Dth_zero(const int k, const int kdth);
        void Dth_minus(const int k, const int kdth);
        void Dth_plus(const int k, const int kdth);
        void Dth_upwind(const int v, const int k, const int kdth);
        void Dthth_zero(const int k, const int kdth);
        
        void Drth_zero(const int k, const int kdrth);
        
        void Dr_zero_4(const int k, const int kdr);
        void Dr_minus_4(const int k, const int kdr);
        void Dr_plus_4(const int k, const int kdr);
        void Dr_upwind_4(const int v, const int kdrm, const int kdrp, const int kdr);
        void Drr_zero_4(const int k, const int kdr);
        
        void Dth_zero_4(const int k, const int kdth);
        void Dth_minus_4(const int k, const int kdth);
        void Dth_plus_4(const int k, const int kdth);
        void Dth_upwind_4(const int v, const int kdrm, const int kdrp, const int kdth);
        void Dthth_zero_4(const int k, const int kdth);
        
        void Drth_zero_4(const int k, const int kdrth);
        
        void NormUpwindGrad(const int v, const int kdxm, const int kdxp,
                            const int kdym, const int kdyp, const int knorm);
        void NormUpwindGrad(const int v, const int k, const int knorm, const int dir=1);
        void NormGrad(const int kdx, const int kdy, const int knorm);
        void Curvature(const int kdx, const int kdy, const int kdxx,
                       const int kdxy, const int kdyy, const int kcurv);
		  
		  void Weight(double w(const double theta), const int kdx, const int kdy,
                    const int kw);

        //single node methods (all one sided difference methods for finding curvature at the particle surface)
		  double Drr_plus(const int i, const int j, const int k);
		  double Drth_plus_zero(const int i, const int j, const int k);
		  double Curvature(const int i, const int j, const double Dr, const double Dt, const double Drr, 
		  						const double Drth, const double Dthth);
		  double Dr_plus_4(const int i, const int j, const int k);		  
		  
		  // Advect methods
        
        double ComputeTVal(const int i, const int j, const int k, const int ikus);
        void UpdateNeighbor(PolarMesh2D_Heap& heap, const int i,
                            const int j, const int k, const int ikus,
                            const int ikhi);

        void ExtendVelocity(const int ktemp, const int kvel, const int ikus,
                                         const int ikhi, const int ikmask);
        //overloaded version allowing initial velocity specified at non-zero phi values
		  void ExtendVelocity(const int kphi, const int ktemp, const int kvel, const int ikus, 
		  												const int ikhi, const int ikmask);
		  void ExtendVals(const int i, const int j, const int k,
                        const int ikus, const int ikhi,
                        const int* p, const int pl);
        void Advance(const int k, const int kv, const int knorm, const double dt);
        void Advance(const int k, const int kv, const double dt);
        
        void Reinitialize(const int k, const int ktemp, const int ik1, const int ik2,
                          const int dir=0);
        void InitLevelSet(const int i, const int j, const int func,
                          const int phi, const int ikus);
    };
}

#endif /* polarmesh2d_h */
