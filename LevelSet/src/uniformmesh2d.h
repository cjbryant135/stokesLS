#ifndef __UNIFORMMESH2D_H__
#define __UNIFORMMESH2D_H__

#ifndef NO_GRAPHICS
#include "plotwindow2d.h"
#include "plotwindow3d.h"
#endif
#include "initialfunc.h"
#include "bicubic.h"
#include "datavec.h"
#include "utility.h"
#include "interface.h"
#include "numtrait.h"
#ifdef USE_AVS
#include "avsfield.h"
#endif
#include <string>
#include <math.h>
#include <limits.h>
#ifdef LEVEL_DEBUG
#include "debug.h"
#endif

// Define boundary conditions for differential operators
// Choices are INWARD (default) and PERIODIC 

// PS - I don't think the above is true...

namespace levelset {
        
#ifndef MESH_ENUMS_DEFINED
#define MESH_ENUMS_DEFINED
    enum STEPSIZE {StepSize};
    enum BOUNDS {Bounds};
    enum Direction {DLeft = 1, DDown = 2, DRight = 4, DUp = 8};
    enum {Upwind = -1, Downwind = 1};
#endif

    class UM2_HeapElement;
    class UniformMesh2D_Heap;
    class UM2_Boundary;

    class UniformMesh2D {
    public:
   
        double      *thedata;            // gridpoint array
        int         datastart;
        int         maxi, maxj, maxk;    // dimensions of data
        int         tmaxi, tmaxj;        // total dimensions of data
        double      dx, dy;              // spatial step sizes
        double      zero[2];             // physical coordinates for i=0,j=0
        int         *theidata;           // integer arrays for grid indexing heap
        int         nextw;
        UM2_Boundary* bc;

        inline int Index(const int i, const int j, const int k=0) const
        {return(k*tmaxi*tmaxj+i*tmaxj+j);}
        inline int Index00(const int i, const int j, const int k=0) const
        {return(k*tmaxi*tmaxj+i*tmaxj+j+datastart);}

    public:

        UniformMesh2D(STEPSIZE, const int m, const int n,
                      const int sd, const int si,
                      const double deltax, const double deltay,
                      const double zx, const double zy,
                      UM2_Boundary& thebc);

        UniformMesh2D(BOUNDS, const int m, const int n,
                      const int sd, const int si,
                      const double xmin, const double xmax,
                      const double ymin, const double ymax,
                      UM2_Boundary& thebc);

        UniformMesh2D(void) {}; 

        virtual ~UniformMesh2D(void);

        // Initialization routines
   
        void SetStepSize(const int m, const int n,
                         const int sd, const int si,
                         const double deltax, const double deltay,
                         const double zx, const double zy, UM2_Boundary& thebc);
        void SetBoundsSize(const int m, const int n,
                           const int sd, const int si,
                           const double xmin, const double xmax, 
                           const double ymin, const double ymax, UM2_Boundary& thebc);
        void SetZero(const double zx, const double zy) {zero[0] = zx; zero[1] = zy;}
   
        void SetValues(const InitialFunc& f, const int wi, const double mult = 1.);
        void InitInterface(const InitialFunc& f, const double start,
                           const double end, const int wi, const int ktemp,
                           const int ik1, const int ik2);
        void SetValues(double (*f)(double,double), const int wi, const double mult = 1.);
        void SetValues(const double* somedata, const int wi, const double mult = 1.);
        void SetValues(const double c, const int wi);
        void InitInterface(double (*fx)(double), double (*fy)(double),
                           const double start, const double end,
                           const int wi);

        // Information retrieval

        inline double X(const int i) const {return zero[0]+i*dx;} 
        inline double Y(const int j) const {return zero[1]+j*dy;} 
        inline double X(const double i) const {return zero[0]+i*dx;} 
        inline double Y(const double j) const {return zero[1]+j*dy;}
        inline int I(const double x) const {return int(floor((x-zero[0])/dx));} 
        inline int J(const double y) const {return int(floor((y-zero[1])/dy));} 
        inline double Id(const double x) const {return (x-zero[0])/dx;} 
        inline double Jd(const double y) const {return (y-zero[1])/dy;} 
        double Interp2(const double& x, const double& y, const int k) const;
        double Interp2(const int i, const int j, const int k, const double& ifrac,
                       const double& jfrac) const;
        Bicubic Interp2(const int i, const int j, const int k) const;

        // Interface only methods
   
        void GetInterface(Interface& segs, const int wi, const double value = 0);
        void UnitNormal(Interface& segs, const int nx, const int ny, 
                        const int kdx, const int kdy);
        void NormUpwindGrad(Interface& segs, const int ng, const int d = 0);
        void Interp(Interface& segs, const int res, const int wi);
        void ExtendVelocity(const Interface& seg, const int sl,
                            const int sh, const int k, const int ktemp,
                            const int ikus, const int ikhi, const int ikmask,
                            const int kvel);
   
		  
		  // Whole mesh methods

        void NormUpwindGrad(const int v, const int kdxm, const int kdxp,
                            const int kdym, const int kdyp, const int knorm);
        void NormUpwindGrad(const int v, const int k, const int knorm, const int dir=1);
        void NormGrad(const int kdx, const int kdy, const int knorm);
        void Curvature(const int kdx, const int kdy, const int kdxx,
                       const int kdxy, const int kdyy, const int kcurv);
        inline void Diffusion(const int kdx, const int kdy, const int kcurv,
                              const int kdiff)
        {Dss_zero(kcurv, kdx, kdy, kdiff);}
        void Weight(double w(const double theta), const int kdx, const int kdy,
                    const int kw);
   
        void Dx_zero(const int k, const int kdx);
        void Dx_minus(const int k, const int kdx);
        void Dx_plus(const int k, const int kdx);
        void Dx_upwind(const int v, const int k, const int kdx);
        void Dxx_zero(const int k, const int kdxx);
        void Dy_zero(const int k, const int kdy);
        void Dy_minus(const int k, const int kdy);
        void Dy_plus(const int k, const int kdy);
        void Dy_upwind(const int v, const int k, const int kdy);
        void Dyy_zero(const int k, const int kdyy);
        void Dxy_zero(const int k, const int kdxy);
        void Dx_zero_4(const int k, const int kdx);
        void Dx_minus_4(const int k, const int kdx);
        void Dx_plus_4(const int k, const int kdx);
        void Dx_upwind_4(const int v, const int kdxm, const int kdxp,
                         const int kdx);
        void Dxx_zero_4(const int k, const int kdxx);
        void Dy_zero_4(const int k, const int kdy);
        void Dy_minus_4(const int k, const int kdy);
        void Dy_plus_4(const int k, const int kdy);
        void Dy_upwind_4(const int v, const int kdym, const int kdyp,
                         const int kdy);
        void Dyy_zero_4(const int k, const int kdyy);
        void Dxy_zero_4(const int k, const int kdxy);
        void Ds_zero(const int kf, const int kdx, const int kdy, const int kds);
        void Ds_zero_4(const int kf, const int kdx, const int kdy, const int kds);
        void Dss_zero(const int kf, const int kdx, const int kdy, const int kds);
        void Dss_zero_4(const int kf, const int kdx, const int kdy, const int kds);

// Single node computations

        double Dx_zero(const int i, const int j, const int k) const ;
        double Dx_minus(const int i, const int j, const int k) const ;
        double Dx_plus(const int i, const int j, const int k) const ;
        double Dx_upwind(const int i, const int j, const int k, const int kv) const ;
        double Dxx_zero(const int i, const int j, const int k) const ;
        double Dy_zero(const int i, const int j, const int k) const ;
        double Dy_minus(const int i, const int j, const int k) const ;
        double Dy_plus(const int i, const int j, const int k) const ;
        double Dy_upwind(const int i, const int j, const int k, const int kv) const ;
        double Dyy_zero(const int i, const int j, const int k) const ;
        double Dxy_zero(const int i, const int j, const int k) const ;
        double Dx_zero_4(const int i, const int j, const int k) const ;
        double Dx_minus_4(const int i, const int j, const int k) const ;
        double Dx_plus_4(const int i, const int j, const int k) const ;
        double Dx_upwind_4(const int i, const int j, const int k, const int kv) const ;
        double Dxx_zero_4(const int i, const int j, const int k) const ;
        double Dy_zero_4(const int i, const int j, const int k) const ;
        double Dy_minus_4(const int i, const int j, const int k) const ;
        double Dy_plus_4(const int i, const int j, const int k) const ;
        double Dy_upwind_4(const int i, const int j, const int k, const int kv) const ;
        double Dyy_zero_4(const int i, const int j, const int k) const ;
        double Dxy_zero_4(const int i, const int j, const int k) const ;
        
        void MinMax(const int k, double& themin, double& themax) const;

        double MaxGrad(const int kdx, const int kdy);
        double AvgGrad(const int kdx, const int kdy);
        double Area(const int k, const double value);
        double IntegrateTrap(const int k);
        double IntegrateBicubic(const int k);
   
        void ExtendVelocity(const int ktemp, const int kvel,
                            const int ikus, const int ikhi, const int ikmask);
        //overloaded for initial velocities at non-zero phi values
        void ExtendVelocity(const int kphi, const int ktemp, const int kvel,
                            const int ikus, const int ikhi, const int ikmask);
		  
		  
		  void Advance(const int k, const int kv, const int knorm, const double dt);
        void Advance(const int k, const int kv, const double dt);
        void Product(const int k, const int ka, const int kb);

        void Reinitialize(const int k, const int ktemp, const int ik1,
                          const int ik2, const int dir = 0);

        char TestValidity(const char* file, const int line,
                          const int kind=-1) const;

        // I/O routines

#ifndef NO_GRAPHICS
        virtual void Plot(PlotWindow2D& port, const int k,
                          const double level_value = 0,
                          const char showboxes = 0x00) const;
        virtual void Plot2(PlotWindow2D& port, const int k,
                           const double level_value = 0, const int reflvl = 10,
                           const char showboxes = 0x00) const;
        virtual void FilledPlot(PlotWindow2D& port, const int k,
                                const unsigned short shade, 
                                const double level_value = 0,
                                const char showboxes = 0x00) const;
        virtual void ShadePlot(PlotWindow2D& Port, const int k, 
                               const double white, const double black, 
                               const char showboxes = 0x00) const;
        virtual void Graph(PlotWindow3D& port, const int k, const double mult = 1.,
                           const char shaded = true) const;
        virtual void ScaledGraph(PlotWindow3D& Port, const int k, const char shaded = true) const;
        virtual void ArrowPlot(PlotWindow2D& port, const int kx, const int ky,
                               const double scale = 1.) const;
        virtual void ArrowPlot(PlotWindow2D& port, const int kx, const int ky,
                               const int kmult, const int kdiv,
                               const double scale = 1.) const;
#endif

        void Plot(std::ostream& s, const int k, const double level_value = 0, 
                  const char showboxes = 0x00) const;
   
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
   
// Utility Methods
   
        char IsCrossing(const double v, const double y0, const double y1,
                        double& f);

        void LocToIndex(const double x, const double y, int& i, int& j,
                        double& ifrac, double& jfrac,
                        const char round = 0x00) const;
        void LocToIndex(const double x, const double y, int& i, int& j, 
                        const char round = 0x00) const;

        inline double Interp(const int i, const int j, const int k, 
                             const double ifrac, const double jfrac) const 
        {return (1.-ifrac)*(1.-jfrac)*data_(i,j,k)
                +(1.-ifrac)*jfrac*data_(i,j+1,k)
                +ifrac*(1.-jfrac)*data_(i+1,j,k)
                +ifrac*jfrac*data_(i+1,j+1,k);}

        inline double Interp(const double& x, const double& y, const int k) const
        {
            int i, j;
            double ifrac, jfrac;
            LocToIndex(x, y, i, j, ifrac, jfrac);
            return Interp(i, j, k, ifrac, jfrac);
        }
   
        void BuildWorkGrid(int num);
        void CopyWorkGrid(int to, int from);
        void CopyIntGrid(int to, int from);

        // Extension Methods
   
        double ComputeTVal(const int i, const int j, const int k);
        double TValFormula(const int sx0, const int sx1, const int sy0, 
                           const int sy1, const double xy0, const double x1, 
                           const double x2, const double y1, const double y2);
        double ComputeTVal(const int i, const int j, const int k, const int ikus);
		  void ExtendDistance(const int phi, const int ikus, const int ikhi, 
                            const int imask);
        void ExtendVals(const int i, const int j, const int k, const int ikus,
                        const int ikhi, const int* p, const int pl);
        void ExtendFunc(const int i, const int j, const int k, const int ikus);
        void UpdateNeighbor(UniformMesh2D_Heap& heap, const int i, const int j, 
                            const int k, const int ikus, const int ikhi);
        void AddInitVelocity(const int i, const int j, const int ikus,
                             const int ikhi, const int func, const int vel,
                             const Direction d, const double f);
        void InitLevelSet(const int i, const int j, const int func, const int phi,
                          const int ikus, const int ikhi);
        void InitLevelSet(const int i, const int j, const int func, const int phi,
                          const int ikus);
        
		  
		  //overloaded velocity extension routines for periodic BCs
        void ExtendVelocity(const int ktemp, const int kvel,
                            const int ikus, const int ikhi, const int ikmask, const bool xper, const bool yper);
        //overloaded for initial velocities at non-zero phi values
        void ExtendVelocity(const int kphi, const int ktemp, const int kvel,
                            const int ikus, const int ikhi, const int ikmask, const bool xper, const bool yper);
		  
		  void UpdateNeighbor(UniformMesh2D_Heap& heap, const int i, const int j, 
                            const int k, const int ikus, const int ikhi, const bool xper, const bool yper, const double phiAcpt);                                
        double ComputeTVal(const int i, const int j, const int k, const int ikus, const bool xper, const bool yper, 
		  							const double phiAcpt);
        void ExtendVals(const int i, const int j, const int k, const int ikus,
                        const int ikhi, const int* p, const int pl, const bool xper, const bool yper);
        
		  // Elliptic equation Methods
   
        void IIM_coeffs(const double x0, const double y0, 
                        const double nx, const double ny, const double chipp, 
                        const double wx, const double wy, const int inside[6], double gamma[6]);
        void Poisson(const int krhs, const int ksoln);

        // I/O methods
   
        friend std::ostream& operator<<(std::ostream& s, UniformMesh2D& m);
        friend std::istream& operator>>(std::istream& s, UniformMesh2D& m);

#ifdef NETCDF
        void ReadNetCDF(const char* fname, const int kstart = 0);
   
        void NetCDF(const char* fname, const int kstart = -1, 
                    const int kend = -1, const char* comment = NULL);
#endif

        void ReadBinary(std::string fname, const int kstart = 0);
        void ReadBinary(const char* fname, const int kstart = 0);
        void WriteBinary(std::string fname, const int kstart = -1,
                         const int kend = -1, const char* comment = NULL);
        void WriteBinary(const char* fname, const int kstart = -1,
                         const int kend = -1, const char* comment = NULL);

#ifdef USE_AVS
        // This function is inline because if not, it doesn't get found by
        // the linker.
        friend ostream& AVS_Out(ostream& s, const UniformMesh2D* m,
                                const int k, const char* comment)
        {
            int dim[2];
            dim[0] = m->maxi;
            dim[1] = m->maxj;
            // NumericTrait<double> t;
            char datatype[20] = "double";
#ifdef AVS_WRITE_UNIFORM
            char field[20] = "uniform";
#else
            char field[20] = "rectilinear";
#endif
            float min_ext[2];
            min_ext[0] = m->zero[0];
            min_ext[1] = m->zero[1];
            float max_ext[2];
            max_ext[0] = m->X(m->maxi-1);
            max_ext[1] = m->Y(m->maxj-1);
            AVSWriteHead(s, comment, 2, dim, 2, 1, datatype, field,
                         min_ext, max_ext, NULL, NULL, NULL, NULL);
            int i, j;
            for (j=0; j<m->maxj; ++j)
                for (i=0; i<m->maxi; ++i)
                    AVSWriteBinary(s, m->mdata_(m,i,j,k));
#ifndef AVS_WRITE_UNIFORM
            for (i=0; i<m->maxi; ++i)
                AVSWriteBinary(s, float(m->X(i)));
            for (j=0; j<m->maxj; ++j)
                AVSWriteBinary(s, float(m->Y(j)));
#endif
            return s;
        }
#endif

        friend class UniformMesh2D_Heap;
        friend class UM2_HeapElement;
   
    };

#define TESTVALIDITY(K) TestValidity(__FILE__,__LINE__,K)

// IO Global operators

    std::ostream& operator<<(std::ostream& s, UniformMesh2D& m);

    std::istream& operator>>(std::istream& s, UniformMesh2D& m);

}

#endif // __UNIFORMMESH2D_H__






