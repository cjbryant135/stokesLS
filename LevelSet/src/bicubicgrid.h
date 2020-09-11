/*************************************************
 bicubicgrid.h
 
 * Revision 1.0  2017/11/28  11:08:23  11:08:23  chopp (David Chopp)
 * Initial revision
 * 
 *************************************************/

#ifndef __BICUBICGRID_H__
#define __BICUBICGRID_H__

#include "uniformmesh2d.h"
#include "bicubic.h"

namespace levelset {

    class BicubicGrid
    {
        const UniformMesh2D* mesh;
        Bicubic** grid;
        int maxi, maxj, kval;
        double zero[2];
        double dx, dy;
        bool periodic[2];
    
    public:
        BicubicGrid(const UniformMesh2D* themesh, const int k, const bool xperiodic, const bool yperiodic);
        ~BicubicGrid(void);

        inline int Index(const int i, const int j) const {return i*maxj+j;}
        inline int I(const double x) const {return floor((x-zero[0])/dx);}
        inline int J(const double y) const {return floor((y-zero[1])/dy);}
        inline double X(const int i) const {return zero[0]+i*dx;}
        inline double Y(const int j) const {return zero[1]+j*dy;}

        double operator()(const double x, const double y);
        double F(const double x, const double y) {return operator()(x,y);}
        double Dx(const double x, const double y);
        double Dy(const double x, const double y);
        double Dxx(const double x, const double y);
        double Dxy(const double x, const double y);
        double Dyy(const double x, const double y);
        double Dxxx(const double x, const double y);
        double Dxxy(const double x, const double y);
        double Dxyy(const double x, const double y);
        double Dyyy(const double x, const double y);
        Bicubic* GetBicubic(const int i, const int j);
        Bicubic* GetBicubic(const double x, const double y);

        double LocalDist(const double x, const double y, 
                         double& ax, double& ay, char& clean);
        double LocalDistNewton(const double x, const double y, 
                               double& ax, double& ay, char& clean);
        double LocalDistBisection(const double x, const double y,
                                  double& ax, double& ay, char& clean);
        double LocalDistDirect(const double x, const double y,
                               double& ax, double& ay);

        void FollowToBdry(const int dir, const double x, const double y, double& ax, double& ay);
    };

}

#endif
