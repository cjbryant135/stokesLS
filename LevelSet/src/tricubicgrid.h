/*************************************************
 tricubicgrid.h
 
 * Revision 1.0  2017/11/28  11:08:23  11:08:23  chopp (David Chopp)
 * Initial revision
 * 
 *************************************************/

#ifndef __TRICUBICGRID_H__
#define __TRICUBICGRID_H__

#include "uniformmesh3d.h"
#include "tricubic.h"

namespace levelset {

    class TricubicGrid
    {
        const UniformMesh3D* mesh;
        Tricubic** grid;
        int maxi, maxj, maxk, kval;
        double zero[3];
        double dx, dy, dz;
        bool periodic[3];
    
    public:
        TricubicGrid(const UniformMesh3D* themesh, const int k, 
                     const bool xperiodic, const bool yperiodic, const bool zperiodic);
        ~TricubicGrid(void);

        inline int Index(const int i, const int j, const int k) const {return (i*maxj+j)*maxk+k;}
        inline int I(const double x) const {return floor((x-zero[0])/dx);}
        inline int J(const double y) const {return floor((y-zero[1])/dy);}
        inline int K(const double z) const {return floor((z-zero[2])/dz);}
        inline double X(const int i) const {return zero[0]+i*dx;}
        inline double Y(const int j) const {return zero[1]+j*dy;}
        inline double Z(const int k) const {return zero[2]+k*dz;}

        double operator()(const double x, const double y, const double z);
        double F(const double x, const double y, const double z) {return operator()(x,y,z);}
        double Dx(const double x, const double y, const double z);
        double Dy(const double x, const double y, const double z);
        double Dz(const double x, const double y, const double z);
        double Dxx(const double x, const double y, const double z);
        double Dxy(const double x, const double y, const double z);
        double Dxz(const double x, const double y, const double z);
        double Dyy(const double x, const double y, const double z);
        double Dyz(const double x, const double y, const double z);
        double Dzz(const double x, const double y, const double z);
        double Dxxx(const double x, const double y, const double z);
        double Dxxy(const double x, const double y, const double z);
        double Dxxz(const double x, const double y, const double z);
        double Dxyy(const double x, const double y, const double z);
        double Dxyz(const double x, const double y, const double z);
        double Dxzz(const double x, const double y, const double z);
        double Dyyy(const double x, const double y, const double z);
        double Dyyz(const double x, const double y, const double z);
        double Dyzz(const double x, const double y, const double z);
        double Dzzz(const double x, const double y, const double z);
        Tricubic* GetTricubic(const int i, const int j, const int k);
        Tricubic* GetTricubic(const double x, const double y, const double z);

        void BestGuess(const double f000, const double f100, const double f010, const double f110,
                       const double f001, const double f101, const double f011, const double f111,
                       const int boxscore, const int corner, double* a) const;
        double LocalDist(const double x, const double y, const double z, 
                         double& ax, double& ay, double& az, char& clean);
        double LocalDistNewton(const double x, const double y, const double z, 
                               double& ax, double& ay, double& az, char& clean);
        double LocalDistBisection(const double x, const double y, const double z,
                                  double& ax, double& ay, double& az, char& clean);

        void FollowToBdry(const int dir, const double x, const double y, const double z, double& ax, double& ay, double& az);
    };

}

#endif
