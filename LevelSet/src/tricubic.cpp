/*************************************************
    tricubic.cpp

    $Header: tricubic.cpp,v 1.1 2000/06/09 22:05:57 chopp Exp $

    $Log:       tricubic.cpp,v $
Revision 1.1  2000/06/09  22:05:57  22:05:57  chopp (David Chopp)
Initial revision

*************************************************/

/*************************************************
    tricubic.cpp

    $Header: tricubic.cpp,v 1.1 2000/06/09 22:05:57 chopp Exp $

    $Log:       tricubic.cpp,v $
Revision 1.1  2000/06/09  22:05:57  22:05:57  chopp (David Chopp)
Initial revision

Revision 1.1  2000/06/01  11:58:41  11:58:41  chopp (David Chopp)
Initial revision

*************************************************/

#include "tricubic.h"
#include <math.h>
#include <iostream>
#include "utility.h"

namespace levelset {
        
    Tricubic::Tricubic(const Tricubic& p)
    {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                for (int k=0; k<4; ++k)
                    a[i][j][k] = p.a[i][j][k];
        wx = p.wx;
        wy = p.wy;
        wz = p.wz;
    }

    void Tricubic::Build(const double f[2][2][2], const double fx[2][2][2], 
                         const double fy[2][2][2], const double fz[2][2][2],
                         const double fxy[2][2][2], const double fxz[2][2][2],
                         const double fyz[2][2][2], const double fxyz[2][2][2],
                         const double dx, const double dy, const double dz)
    {
        int i, j, k;
        double g[64];
        wx = dx;
        wy = dy;
        wz = dz;

        for (i=0; i<=1; ++i)
            for (j=0; j<=1; ++j)
                for ( k=0; k<=1; ++k) {
                    g[4*i+2*j+k] = f[i][j][k];
                    g[4*i+2*j+k+8] = fx[i][j][k]*dx;
                    g[4*i+2*j+k+16] = fy[i][j][k]*dy;
                    g[4*i+2*j+k+24] = fz[i][j][k]*dz;
                    g[4*i+2*j+k+32] = fxy[i][j][k]*dx*dy;
                    g[4*i+2*j+k+40] = fxz[i][j][k]*dx*dz;
                    g[4*i+2*j+k+48] = fyz[i][j][k]*dy*dz;
                    g[4*i+2*j+k+56] = fxyz[i][j][k]*dx*dy*dz;
                }   

        static int A[64][64] = {
            1, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  1, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,
      
            -3, 3,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0, -2,-1,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,
      
            2,-2,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  1, 1,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,
      
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  1, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  1, 0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,
      
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -3, 3, 0,0,  0, 0, 0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0, -2,-1, 0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,
      
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  2,-2, 0,0,  0, 0, 0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  1, 1, 0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            -3, 0,3,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -2, 0,-1,0,  0,0, 0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0,0, 0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0,  0,0,0,0, -3, 0,3,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0, 0,0, -2,0,-1,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            9,-9,-9,9, 0,0,0,0, 0,0,0,0, 0,0,0,0, 6,-6,3,-3, 0,0,0,0,  6,3,-6,-3, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 4,2,2,1,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            -6,6,6,-6, 0,0,0,0, 0,0,0,0, 0,0,0,0, -4,4,-2,2, 0,0,0,0,  -3,-3,3,3, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -2,-2,-1,-1,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            2, 0,-2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,1,0, 0,0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  2,0,-2,0, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,1,0,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            -6,6,6,-6, 0,0,0,0, 0,0,0,0, 0,0,0,0, -3,3,-3,3, 0,0,0,0,  -4,-2,4,2, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -2,-1,-2,-1,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            4,-4,-4,4, 0,0,0,0, 0,0,0,0, 0,0,0,0, 2,-2,2,-2, 0,0,0,0,  2,2,-2,-2, 
            0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1,  0, 0,0,0, 
            0, 0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, -3,3,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, -2,-1,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 2,-2,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            1,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, -3,3,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            -2,-1,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 2,-2,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            1,1,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, -3,0,3,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, -2,0,-1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, -3,0,3,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            -2,0,-1,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 9,-9,-9,9, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 6,-6,3,-3, 0,0,0,0, 6,3,-6,-3, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            4,2,2,1, 0,0,0,0,

            0,0,0,0, 0,0,0,0, -6,6,6,-6, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, -4,4,-2,2, 0,0,0,0, -3,-3,3,3, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            -2,-2,-1,-1, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 2,0,-2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 1,0,1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 2,0,-2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            1,0,1,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, -6,6,6,-6, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, -3,3,-3,3, 0,0,0,0, -4,-2,4,2, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            -2,-1,-2,-1, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 4,-4,-4,4, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 2,-2,2,-2, 0,0,0,0, 2,2,-2,-2, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            1,1,1,1, 0,0,0,0,

            -3,0,0,0, 3,0,0,0, -2,0,0,0, -1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -3,0,0,0, 
            3,0,0,0, 0,0,0,0, 0,0,0,0, -2,0,0,0, -1,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            9,-9,0,0, -9,9,0,0, 6,-6,0,0, 3,-3,0,0, 0,0,0,0, 0,0,0,0, 6,3,0,0, 
            -6,-3,0,0, 0,0,0,0, 0,0,0,0, 4,2,0,0, 2,1,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            -6,6,0,0, 6,-6,0,0, -4,4,0,0, -2,2,0,0, 0,0,0,0, 0,0,0,0, -3,-3,0,0, 
            3,3,0,0, 0,0,0,0, 0,0,0,0, -2,-2,0,0, -1,-1,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -3,0,0,0, 3,0,0,0, 0,0,0,0, 
            0,0,0,0, -2,0,0,0, -1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -3,0,0,0, 3,0,0,0, 
            -2,0,0,0, -1,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 9,-9,0,0, -9,9,0,0, 0,0,0,0, 
            0,0,0,0, 6,-6,0,0, 3,-3,0,0, 0,0,0,0, 0,0,0,0, 6,3,0,0, -6,-3,0,0, 
            4,2,0,0, 2,1,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -6,6,0,0, 6,-6,0,0, 0,0,0,0, 
            0,0,0,0, -4,4,0,0, -2,2,0,0, 0,0,0,0, 0,0,0,0, -3,-3,0,0, 3,3,0,0, 
            -2,-2,0,0, -1,-1,0,0,

            9,0,-9,0, -9,0,9,0, 6,0,-6,0, 3,0,-3,0, 6,0,3,0, -6,0,-3,0, 0,0,0,0, 
            0,0,0,0, 4,0,2,0, 2,0,1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 9,0,-9,0, 
            -9,0,9,0, 0,0,0,0, 0,0,0,0, 6,0,-6,0, 3,0,-3,0, 6,0,3,0, -6,0,-3,0, 
            4,0,2,0, 2,0,1,0,

            -27,27,27,-27, 27,-27,-27,27, -18,18,18,-18, -9,9,9,-9, -18,18,-9,9, 
            18,-18,9,-9, -18,-9,18,9, 18,9,-18,-9, -12,12,-6,6, -6,6,-3,3,
            -12,-6,12,6, -6,-3,6,3, -12,-6,-6,-3, 12,6,6,3, -8,-4,-4,-2,
            -4,-2,-2,-1,

            18,-18,-18,18, -18,18,18,-18, 12,-12,-12,12, 6,-6,-6,6, 12,-12,6,-6,
            -12,12,-6,6, 9,9,-9,-9, -9,-9,9,9, 8,-8,4,-4, 4,-4,2,-2, 6,6,-6,-6,
            3,3,-3,-3, 6,6,3,3, -6,-6,-3,-3, 4,4,2,2, 2,2,1,1,

            -6,0,6,0, 6,0,-6,0, -4,0,4,0, -2,0,2,0, -3,0,-3,0, 3,0,3,0, 0,0,0,0, 
            0,0,0,0, -2,0,-2,0, -1,0,-1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -6,0,6,0, 
            6,0,-6,0, 0,0,0,0, 0,0,0,0, -4,0,4,0, -2,0,2,0, -3,0,-3,0, 3,0,3,0, 
            -2,0,-2,0, -1,0,-1,0,

            18,-18,-18,18, -18,18,18,-18, 12,-12,-12,12, 6,-6,-6,6, 9,-9,9,-9,
            -9,9,-9,9, 12,6,-12,-6, -12,-6,12,6, 6,-6,6,-6, 3,-3,3,-3, 8,4,-8,-4,
            4,2,-4,-2, 6,3,6,3, -6,-3,-6,-3, 4,2,4,2, 2,1,2,1, 

            -12,12,12,-12, 12,-12,-12,12, -8,8,8,-8, -4,4,4,-4, -6,6,-6,6,
            6,-6,6,-6, -6,-6,6,6, 6,6,-6,-6, -4,4,-4,4, -2,2,-2,2, -4,-4,4,4,
            -2,-2,2,2, -3,-3,-3,-3, 3,3,3,3, -2,-2,-2,-2, -1,-1,-1,-1,

            2,0,0,0, -2,0,0,0, 1,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 2,0,0,0, 
            -2,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            -6,6,0,0, 6,-6,0,0, -3,3,0,0, -3,3,0,0, 0,0,0,0, 0,0,0,0, -4,-2,0,0, 
            4,2,0,0, 0,0,0,0, 0,0,0,0, -2,-1,0,0, -2,-1,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            4,-4,0,0, -4,4,0,0, 2,-2,0,0, 2,-2,0,0, 0,0,0,0, 0,0,0,0, 2,2,0,0, 
            -2,-2,0,0, 0,0,0,0, 0,0,0,0, 1,1,0,0, 1,1,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 2,0,0,0, -2,0,0,0, 0,0,0,0, 
            0,0,0,0, 1,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 2,0,0,0, -2,0,0,0, 
            1,0,0,0, 1,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -6,6,0,0, 6,-6,0,0, 0,0,0,0, 
            0,0,0,0, -3,3,0,0, -3,3,0,0, 0,0,0,0, 0,0,0,0, -4,-2,0,0, 4,2,0,0, 
            -2,-1,0,0, -2,-1,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 4,-4,0,0, -4,4,0,0, 0,0,0,0, 
            0,0,0,0, 2,-2,0,0, 2,-2,0,0, 0,0,0,0, 0,0,0,0, 2,2,0,0, -2,-2,0,0, 
            1,1,0,0, 1,1,0,0,

            -6,0,6,0, 6,0,-6,0, -3,0,3,0, -3,0,3,0, -4,0,-2,0, 4,0,2,0, 0,0,0,0, 
            0,0,0,0, -2,0,-1,0, -2,0,-1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, -6,0,6,0, 
            6,0,-6,0, 0,0,0,0, 0,0,0,0, -3,0,3,0, -3,0,3,0, -4,0,-2,0, 4,0,2,0, 
            -2,0,-1,0, -2,0,-1,0,

            18,-18,-18,18, -18,18,18,-18, 9,-9,-9,9, 9,-9,-9,9, 12,-12,6,-6,
            -12,12,-6,6, 12,6,-12,-6, -12,-6,12,6, 6,-6,3,-3, 6,-6,3,-3,
            6,3,-6,-3, 6,3,-6,-3, 8,4,4,2, -8,-4,-4,-2, 4,2,2,1, 4,2,2,1,

            -12,12,12,-12, 12,-12,-12,12, -6,6,6,-6, -6,6,6,-6, -8,8,-4,4,
            8,-8,4,-4, -6,-6,6,6, 6,6,-6,-6, -4,4,-2,2, -4,4,-2,2, -3,-3,3,3,
            -3,-3,3,3, -4,-4,-2,-2, 4,4,2,2, -2,-2,-1,-1, -2,-2,-1,-1,

            4,0,-4,0, -4,0,4,0, 2,0,-2,0, 2,0,-2,0, 2,0,2,0, -2,0,-2,0, 0,0,0,0, 
            0,0,0,0, 1,0,1,0, 1,0,1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
            0,0,0,0, 0,0,0,0,

            0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 4,0,-4,0, 
            -4,0,4,0, 0,0,0,0, 0,0,0,0, 2,0,-2,0, 2,0,-2,0, 2,0,2,0, -2,0,-2,0, 
            1,0,1,0, 1,0,1,0,

            -12,12,12,-12, 12,-12,-12,12, -6,6,6,-6, -6,6,6,-6, -6,6,-6,6,
            6,-6,6,-6, -8,-4,8,4, 8,4,-8,-4, -3,3,-3,3, -3,3,-3,3, -4,-2,4,2,
            -4,-2,4,2, -4,-2,-4,-2, 4,2,4,2, -2,-1,-2,-1, -2,-1,-2,-1,

            8,-8,-8,8, -8,8,8,-8, 4,-4,-4,4, 4,-4,-4,4, 4,-4,4,-4, -4,4,-4,4,
            4,4,-4,-4, -4,-4,4,4, 2,-2,2,-2, 2,-2,2,-2, 2,2,-2,-2, 2,2,-2,-2,
            2,2,2,2, -2,-2,-2,-2, 1,1,1,1, 1,1,1,1
        };

        int l,n;
   
        for (i=0, n=0; i<4; ++i)
            for (j=0; j<4; ++j)
                for (k=0; k<4; ++k, ++n) {
                    a[i][j][k] = 0.;
                    double denom = pow(wx,i)*pow(wy,j)*pow(wz,k); 
                    for (l=0; l<64; ++l)
                        a[i][j][k] += A[n][l]*g[l]/denom;
                }   
    }

    void Tricubic::BuildwDeriv(const double f[4][4][4], const double dx,
                               const double dy, const double dz)
    {
        double fu[2][2][2];
        double fx[2][2][2];
        double fy[2][2][2];
        double fz[2][2][2];
        double fxy[2][2][2];
        double fxz[2][2][2];
        double fyz[2][2][2];
        double fxyz[2][2][2];
        
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j)
                for (int k=0; k<2; ++k) {
                    fu[i][j][k] = f[1+i][1+j][1+k];
                    fx[i][j][k] = (f[2+i][1+j][1+k]-f[i][1+j][1+k])/2./dx;
                    fy[i][j][k] = (f[1+i][2+j][1+k]-f[1+i][j][1+k])/2./dy;
                    fz[i][j][k] = (f[1+i][1+j][2+k]-f[1+i][1+j][k])/2./dz;
                    fxy[i][j][k] = (f[2+i][2+j][1+k]-f[i][2+j][1+k]
                                    -f[2+i][j][1+k]+f[i][j][1+k])/4./dx/dy;
                    fxz[i][j][k] = (f[2+i][1+j][2+k]-f[i][1+j][2+k]
                                    -f[2+i][1+j][k]+f[i][1+j][k])/4./dx/dz;
                    fyz[i][j][k] = (f[1+i][2+j][2+k]-f[1+i][j][2+k]
                                    -f[1+i][2+j][k]+f[1+i][j][k])/4./dy/dz;
                    fxyz[i][j][k] = (f[2+i][2+j][2+k]-f[i][2+j][2+k]-f[2+i][j][2+k]+f[i][j][2+k]
                                     -f[2+i][2+j][k]+f[i][2+j][k]+f[2+i][j][k]-f[i][j][k])/8./dx/dy/dz;
                }
                        
        Build(fu, fx, fy, fz, fxy, fxz, fyz, fxyz, dx, dy, dz);
    }

    double Tricubic::operator()(const double x, const double y, 
                                const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += a[i][j][k]*pow(x,i)*pow(y,j)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dx(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=1; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += i*a[i][j][k]*pow(x,i-1)*pow(y,j)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dy(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=1; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += j*a[i][j][k]*pow(x,i)*pow(y,j-1)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=1; k<4; ++k)
                    ans += k*a[i][j][k]*pow(x,i)*pow(y,j)*pow(z,k-1);
    
        return ans;
    }

    double Tricubic::Dxx(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=2; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += i*(i-1)*a[i][j][k]*pow(x,i-2)*pow(y,j)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dxy(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=1; i<4; ++i) 
            for (int j=1; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += i*j*a[i][j][k]*pow(x,i-1)*pow(y,j-1)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dxz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=1; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=1; k<4; ++k)
                    ans += i*k*a[i][j][k]*pow(x,i-1)*pow(y,j)*pow(z,k-1);
    
        return ans;
    }

    double Tricubic::Dyy(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=2; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += j*(j-1)*a[i][j][k]*pow(x,i)*pow(y,j-2)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dyz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=1; j<4; ++j)
                for (int k=1; k<4; ++k)
                    ans += j*k*a[i][j][k]*pow(x,i)*pow(y,j-1)*pow(z,k-1);
    
        return ans;
    }

    double Tricubic::Dzz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=2; k<4; ++k)
                    ans += k*(k-1)*a[i][j][k]*pow(x,i)*pow(y,j)*pow(z,k-2);
    
        return ans;
    }

    double Tricubic::Dxxx(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int j=0; j<4; ++j)
            for (int k=0; k<4; ++k)
                ans += 6*a[3][j][k]*pow(y,j)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dxxy(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=2; i<4; ++i) 
            for (int j=1; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += i*(i-1)*j*a[i][j][k]*pow(x,i-2)*pow(y,j-1)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dxxz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=2; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=1; k<4; ++k)
                    ans += i*(i-1)*k*a[i][j][k]*pow(x,i-1)*pow(y,j)*pow(z,k-1);
    
        return ans;
    }

    double Tricubic::Dxyy(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=1; i<4; ++i) 
            for (int j=2; j<4; ++j)
                for (int k=0; k<4; ++k)
                    ans += i*j*(j-1)*a[i][j][k]*pow(x,i-1)*pow(y,j-2)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dxyz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=1; i<4; ++i) 
            for (int j=1; j<4; ++j)
                for (int k=1; k<4; ++k)
                    ans += i*j*k*a[i][j][k]*pow(x,i-1)*pow(y,j-1)*pow(z,k-1);
    
        return ans;
    }

    double Tricubic::Dxzz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=1; i<4; ++i) 
            for (int j=0; j<4; ++j)
                for (int k=2; k<4; ++k)
                    ans += i*k*(k-1)*a[i][j][k]*pow(x,i-1)*pow(y,j)*pow(z,k-2);
    
        return ans;
    }

    double Tricubic::Dyyy(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int k=0; k<4; ++k)
                ans += 6*a[i][3][k]*pow(x,i)*pow(z,k);
    
        return ans;
    }

    double Tricubic::Dyyz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=2; j<4; ++j)
                for (int k=1; k<4; ++k)
                    ans += j*(j-1)*k*a[i][j][k]*pow(x,i)*pow(y,j-2)*pow(z,k-1);
    
        return ans;
    }

    double Tricubic::Dyzz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=1; j<4; ++j)
                for (int k=2; k<4; ++k)
                    ans += j*k*(k-1)*a[i][j][k]*pow(x,i)*pow(y,j-1)*pow(z,k-2);
    
        return ans;
    }

    double Tricubic::Dzzz(const double x, const double y, const double z) const
    {
        double ans = 0;
    
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j)
                ans += 6*a[i][j][3]*pow(x,i)*pow(y,j);
    
        return ans;
    }

    inline double sqr(const double x) {return x*x;}

#define MAXIT 40

    double Tricubic::LocalDist(const double x, const double y, 
                               const double z, double& ax, double& ay, 
                               double& az, char& clean, const double tol) const
    {
        ax = x;
        ay = y;
        az = z;
        double p, px, py, pz, tax, tay, taz;
        double delta[6] = {1000.,1000.,1000.,1000.,1000.,1000.};
        int itcount = 0;

        while (sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3])
                    +sqr(delta[4])+sqr(delta[5])) > tol*wx*wy 
               && itcount < MAXIT) {

            // Find zero

            p = F(ax,ay,az);
            px = Dx(ax,ay,az);
            py = Dy(ax,ay,az);
            pz = Dz(ax,ay,az);
            delta[0] = -p*px/(px*px+py*py+pz*pz);
            delta[1] = -p*py/(px*px+py*py+pz*pz);
            delta[2] = -p*pz/(px*px+py*py+pz*pz);
            tax = ax+delta[0];
            tay = ay+delta[1];
            taz = az+delta[2];

            // Find min distance

            p = F(tax,tay,taz);
            px = Dx(tax,tay,taz);
            py = Dy(tax,tay,taz);
            pz = Dz(tax,tay,taz);
            delta[3] = ((py*py+pz*pz)*(x-tax)-px*(py*(y-tay)+pz*(z-taz)))
                /(px*px+py*py+pz*pz);
            delta[4] = ((px*px+pz*pz)*(y-tay)-py*(px*(x-tax)+pz*(z-taz)))
                /(px*px+py*py+pz*pz);
            delta[5] = ((px*px+py*py)*(z-taz)-pz*(px*(x-tax)+py*(y-tay)))
                /(px*px+py*py+pz*pz);
      
            ax += delta[0]+delta[3];
            ay += delta[1]+delta[4];
            az += delta[2]+delta[5];
            ++itcount;
        }

#if 0
        cout << "Iterations = " << itcount << '\n';
#endif

        clean = itcount < MAXIT && 0 <= ax && ax <= wx && 0 <= ay && ay <= wy
            && 0 <= az && az <= wz;
   
        return sqrt(sqr(x-ax)+sqr(y-ay)+sqr(z-az));
    }

    Bicubic Tricubic::Slice(const double xyz, const int dir) const
    {
        Bicubic p;
        
        switch(dir) {
        case 0:    
            for (int j=0; j<4; ++j)
                for (int k=0; k<4; ++k) {
                    p.a[j][k] = 0.;
                    for (int i=0; i<4; ++i) 
                        p.a[j][k] += a[i][j][k]*pow(xyz,i);
                }
            p.wx = wy;
            p.wy = wz;
            break;
        case 1:
            for (int i=0; i<4; ++i)
                for (int k=0; k<4; ++k) {
                    p.a[i][k] = 0.;
                    for (int j=0; j<4; ++j)
                        p.a[i][k] += a[i][j][k]*pow(xyz,j);
                }
            p.wx = wx;
            p.wy = wz;
            break;
        case 2:
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    p.a[i][j] = 0.;
                    for (int k=0; k<4; ++k)
                        p.a[i][j] += a[i][j][k]*pow(xyz,k);
                }
            p.wx = wx;
            p.wy = wy;
            break;
        }
        return p;
    }


    void Tricubic::XCrossing(const double x, const double y, const double z, 
                             double& ay, double& az, char& clean) const 
    {
        double yy, zz;
        Bicubic p = Slice(x, 0);
        p.LocalDist(y, z, ay, az, clean);
        if (ay < 0) {
            if (az < 0) {
                zz = p.XCrossing(0., z, zz, clean);
                yy = p.YCrossing(y, 0., yy, clean);
                if (sqrt(sqr(y-yy)+z*z) < sqrt(y*y+sqr(z-zz))) {
                    ay = yy;
                    az = 0.;
                } else {
                    ay = 0.;
                    az = zz;
                }
            } else if (az > wz) {
                zz = p.XCrossing(0., z, zz, clean);
                yy = p.YCrossing(y, wz, yy, clean);
                if (sqrt(sqr(y-yy)+sqr(z-wz)) < sqrt(y*y+sqr(z-zz))) {
                    ay = yy;
                    az = wz;
                } else {
                    ay = 0.;
                    az = zz;
                }
            } else {
                ay = 0.;
                zz = p.XCrossing(0., z, zz, clean);
            }
        } else if (ay > wy) {
            if (az < 0) {
                zz = p.XCrossing(wy, z, zz, clean);
                yy = p.YCrossing(y, 0., zz, clean);
                if (sqrt(sqr(y-yy)+sqr(z)) < sqrt(sqr(y)+sqr(z-zz))) {
                    ay = yy;
                    az = 0.;
                } else {
                    ay = wy;
                    az = zz;
                }
            } else if (az > wz) {
                zz = p.XCrossing(wy, z, zz, clean);
                yy = p.YCrossing(y, wz, yy, clean);
                if (sqrt(sqr(y-yy)+sqr(z-wz)) < sqrt(sqr(y-wy)+sqr(z-zz))) {
                    ay = yy;
                    az = wz;
                } else {
                    ay = wy;
                    az = zz;
                }
            } else {
                ay = wy;
                zz = p.XCrossing(wy, z, zz, clean);
            }
        } else {
            if (az < 0) {
                ay = p.YCrossing(y, 0., yy, clean);
                az = 0.;
            } else if (az > wz) {
                ay = p.YCrossing(y, wz, yy, clean);
                az = wz;
            }
        }
    }

    void Tricubic::YCrossing(const double x, const double y, const double z, 
                             double& ax, double& az, char& clean) const 
    {
        double xx, zz;
        Bicubic p = Slice(y, 1);
        p.LocalDist(x, z, ax, az, clean);
        if (ax < 0) {
            if (az < 0) {
                zz = p.XCrossing(0., z, zz, clean);
                xx = p.YCrossing(x, 0., xx, clean);
                if (sqrt(sqr(x-xx)+z*z) < sqrt(x*x+sqr(z-zz))) {
                    ax = xx;
                    az = 0.;
                } else {
                    ax = 0.;
                    az = zz;
                }
            } else if (az > wz) {
                zz = p.XCrossing(0., z, zz, clean);
                xx = p.YCrossing(x, wz, xx, clean);
                if (sqrt(sqr(x-xx)+sqr(z-wz)) < sqrt(x*x+sqr(z-zz))) {
                    ax = xx;
                    az = wz;
                } else {
                    ax = 0.;
                    az = zz;
                }
            } else {
                ax = 0.;
                zz = p.XCrossing(0., z, zz, clean);
            }
        } else if (ax > wx) {
            if (az < 0) {
                zz = p.XCrossing(wx, z, zz, clean);
                xx = p.YCrossing(x, 0., zz, clean);
                if (sqrt(sqr(x-xx)+sqr(z)) < sqrt(sqr(x)+sqr(z-zz))) {
                    ax = xx;
                    az = 0.;
                } else {
                    ax = wx;
                    az = zz;
                }
            } else if (az > wz) {
                zz = p.XCrossing(wx, z, zz, clean);
                xx = p.YCrossing(x, wz, xx, clean);
                if (sqrt(sqr(x-xx)+sqr(z-wz)) < sqrt(sqr(x-wx)+sqr(z-zz))) {
                    ax = xx;
                    az = wz;
                } else {
                    ax = wx;
                    az = zz;
                }
            } else {
                ax = wx;
                zz = p.XCrossing(wx, z, zz, clean);
            }
        } else {
            if (az < 0) {
                ax = p.YCrossing(x, 0., xx, clean);
                az = 0.;
            } else if (az > wz) {
                ax = p.YCrossing(x, wz, xx, clean);
                az = wz;
            }
        }
    }

    void Tricubic::ZCrossing(const double x, const double y, const double z, 
                             double& ax, double& ay, char& clean) const
    {
        double xx, yy;
        Bicubic p = Slice(y, 1);
        p.LocalDist(x, y, ax, ay, clean);
        if (ax < 0) {
            if (ay < 0) {
                yy = p.XCrossing(0., y, yy, clean);
                xx = p.YCrossing(x, 0., xx, clean);
                if (sqrt(sqr(x-xx)+y*y) < sqrt(x*x+sqr(y-yy))) {
                    ax = xx;
                    ay = 0.;
                } else {
                    ax = 0.;
                    ay = yy;
                }
            } else if (ay > wy) {
                yy = p.XCrossing(0., y, yy, clean);
                xx = p.YCrossing(x, wy, xx, clean);
                if (sqrt(sqr(x-xx)+sqr(y-wy)) < sqrt(x*x+sqr(y-yy))) {
                    ax = xx;
                    ay = wy;
                } else {
                    ax = 0.;
                    ay = yy;
                }
            } else {
                ax = 0.;
                yy = p.XCrossing(0., y, yy, clean);
            }
        } else if (ax > wx) {
            if (ay < 0) {
                yy = p.XCrossing(wx, y, yy, clean);
                xx = p.YCrossing(x, 0., yy, clean);
                if (sqrt(sqr(x-xx)+sqr(y)) < sqrt(sqr(x)+sqr(y-yy))) {
                    ax = xx;
                    ay = 0.;
                } else {
                    ax = wx;
                    ay = yy;
                }
            } else if (ay > wy) {
                yy = p.XCrossing(wx, y, yy, clean);
                xx = p.YCrossing(x, wy, xx, clean);
                if (sqrt(sqr(x-xx)+sqr(y-wy)) < sqrt(sqr(x-wx)+sqr(y-yy))) {
                    ax = xx;
                    ay = wy;
                } else {
                    ax = wx;
                    ay = yy;
                }
            } else {
                ax = wx;
                yy = p.XCrossing(wx, y, yy, clean);
            }
        } else {
            if (ay < 0) {
                ax = p.YCrossing(x, 0., xx, clean);
                ay = 0.;
            } else if (ay > wy) {
                ax = p.YCrossing(x, wy, xx, clean);
                ay = wy;
            }
        }
    }

}


