//
//  genbicubic.cpp
//  LevelSet
//
//  Created by David Chopp on 8/8/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "genbicubic.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "lapack.h"

namespace levelset {
    
    //    inline double mypow(const double a, const double b) {return a != 0. ? pow(a,b) : 0.;}
    inline double mypow(const double a, const double b) {return pow(a,b);}
    
    GeneralBicubic::GeneralBicubic(const GeneralBicubic& p)
    {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                a[i][j] = p.a[i][j];
        for (int i=0; i<3; ++i) {
            wx[i] = p.wx[i];
            wy[i] = p.wy[i];
        }
    }
    
    
    void GeneralBicubic::Build(const double f11, const double f21, const double f12,
                               const double f22, const double fx11, const double fx21,
                               const double fx12, const double fx22, const double fy11,
                               const double fy21, const double fy12, const double fy22,
                               const double fxy11, const double fxy21, const double fxy12,
                               const double fxy22, const double dx0, const double dx1,
                               const double dx2, const double dy0, const double dy1,
                               const double dy2)
    {
        double f[16];
        f[0] = f11;
        f[1] = f21;
        f[2] = f12;
        f[3] = f22;
        f[4] = fx11*dx1;
        f[5] = fx21*dx1;
        f[6] = fx12*dx1;
        f[7] = fx22*dx1;
        f[8] = fy11*dy1;
        f[9] = fy21*dy1;
        f[10] = fy12*dy1;
        f[11] = fy22*dy1;
        f[12] = fxy11*dx1*dy1;
        f[13] = fxy21*dx1*dy1;
        f[14] = fxy12*dx1*dy1;
        f[15] = fxy22*dx1*dy1;
        
        a[0][0] = f[0];
        a[0][1] = f[8]/dy1;
        a[0][2] = (-3*f[0]+3*f[2]-2*f[8]-f[10])/dy1/dy1;
        a[0][3] = (2*f[0]-2*f[2]+f[8]+f[10])/dy1/dy1/dy1;
        a[1][0] = f[4]/dx1;
        a[1][1] = f[12]/dx1/dy1;
        a[1][2] = (-3*f[4]+3*f[6]-2*f[12]-f[14])/dx1/dy1/dy1;
        a[1][3] = (2*f[4]-2*f[6]+f[12]+f[14])/dx1/dy1/dy1/dy1;
        a[2][0] = (-3*f[0]+3*f[1]-2*f[4]-f[5])/dx1/dx1;
        a[2][1] = (-3*f[8]+3*f[9]-2*f[12]-f[13])/dx1/dx1/dy1;
        a[2][2] =  (9*f[0]-9*f[1]-9*f[2]+9*f[3]+6*f[4]+3*f[5]-6*f[6]-3*f[7]+6*f[8]
                    -6*f[9]+3*f[10]-3*f[11]+4*f[12]+2*f[13]+2*f[14]+f[15])
        /dx1/dx1/dy1/dy1;
        a[2][3] = (-6*f[0]+6*f[1]+6*f[2]-6*f[3]-4*f[4]-2*f[5]+4*f[6]+2*f[7]-3*f[8]
                   +3*f[9]-3*f[10]+3*f[11]-2*f[12]-f[13]-2*f[14]-f[15])
        /dx1/dx1/dy1/dy1/dy1;
        a[3][0] =  (2*f[0]-2*f[1]+f[4]+f[5])/dx1/dx1/dx1;
        a[3][1] =  (2*f[8]-2*f[9]+f[12]+f[13])/dx1/dx1/dx1/dy1;
        a[3][2] = (-6*f[0]+6*f[1]+6*f[2]-6*f[3]-3*f[4]-3*f[5]+3*f[6]+3*f[7]-4*f[8]
                   +4*f[9]-2*f[10]+2*f[11]-2*f[12]-2*f[13]-f[14]-f[15])
        /dx1/dx1/dx1/dy1/dy1;
        a[3][3] =  (4*f[0]-4*f[1]-4*f[2]+4*f[3]+2*f[4]+2*f[5]-2*f[6]-2*f[7]+2*f[8]
                    -2*f[9]+2*f[10]-2*f[11]+f[12]+f[13]+f[14]+f[15])
        /dx1/dx1/dx1/dy1/dy1/dy1;
        wx[0] = dx0; wx[1] = dx1; wx[2] = dx2;
        wy[0] = dy0; wy[1] = dy1; wy[2] = dy2;
    }
    
    void GeneralBicubic::BuildwDeriv(const double f00, const double f10, const double f20,
                                     const double f30, const double f01, const double f11,
                                     const double f21, const double f31, const double f02,
                                     const double f12, const double f22, const double f32,
                                     const double f03, const double f13, const double f23,
                                     const double f33, const double dx0, const double dx1,
                                     const double dx2, const double dy0, const double dy1,
                                     const double dy2)
    {
#if 0
        const double A[16][16] = {{  0.,   0.,   0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
            {  0.,   0.,   0.,   0., -0.5,   0.,  0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
            {  0.,   0.,   0.,   0.,   1., -2.5,   2., -0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
            {  0.,   0.,   0.,   0., -0.5,  1.5, -1.5,  0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
            {  0., -0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5,   0.,   0.,   0.,   0.,   0.,   0.},
            {0.25,   0.,-0.25,   0.,   0.,   0.,   0.,   0.,-0.25,   0., 0.25,   0.,   0.,   0.,   0.,   0.},
            {-0.5, 1.25,  -1., 0.25,   0.,   0.,   0.,   0.,  0.5,-1.25,   1.,-0.25,   0.,   0.,   0.,   0.},
            {0.25,-0.75, 0.75,-0.25,   0.,   0.,   0.,   0.,-0.25, 0.75,-0.75, 0.25,   0.,   0.,   0.,   0.},
            {  0.,   1.,   0.,   0.,   0., -2.5,   0.,   0.,   0.,   2.,   0.,   0.,   0., -0.5,   0.,   0.},
            {-0.5,   0.,  0.5,   0., 1.25,   0.,-1.25,   0.,  -1.,   0.,   1.,   0., 0.25,   0.,-0.25,   0.},
            {  1., -2.5,   2., -0.5, -2.5, 6.25,  -5., 1.25,   2.,  -5.,   4.,  -1., -0.5, 1.25,  -1., 0.25},
            {-0.5,  1.5, -1.5,  0.5, 1.25,-3.75, 3.75,-1.25,  -1.,   3.,  -3.,   1., 0.25,-0.75, 0.75,-0.25},
            {  0., -0.5,   0.,   0.,   0.,  1.5,   0.,   0.,   0., -1.5,   0.,   0.,   0.,  0.5,   0.,   0.},
            {0.25,   0.,-0.25,   0.,-0.75,   0., 0.75,   0., 0.75,   0.,-0.75,   0.,-0.25,   0., 0.25,   0.},
            {-0.5, 1.25,  -1., 0.25,  1.5,-3.75,   3.,-0.75, -1.5, 3.75,  -3., 0.75,  0.5,-1.25,   1.,-0.25},
            {0.25,-0.75, 0.75,-0.25,-0.75, 2.25,-2.25, 0.75, 0.75,-2.25, 2.25,-0.75,-0.25, 0.75,-0.75, 0.25}};
#endif
        
        const double A[16][16] = {
            {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, -(dy1/(dy0*(dy0 + dy1))), 1/dy0 - 1/dy1,
                dy0/(dy1*(dy0 + dy1)), 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 2/(dy0*(dy0 + dy1)), -((2/dy0 + 1/(dy1 + dy2))/dy1),
                (dy0 + dy1 + 2*dy2)/(dy0*dy1*dy2 + dy1*dy1*dy2),
                -(1/(dy2*(dy1 + dy2))), 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, -(1/(dy0*dy1*(dy0 + dy1))),
                (dy0 + dy1 + dy2)/(dy0*dy1*dy1*(dy1 + dy2)),
                -((dy0 + dy1 + dy2)/(dy1*dy1*(dy0 + dy1)*dy2)),
                1/(dy1*dy1*dy2 + dy1*dy2*dy2), 0, 0, 0, 0, 0, 0, 0, 0},
            {0, -(dx1/(dx0*(dx0 + dx1))), 0, 0, 0, 1/dx0 - 1/dx1, 0, 0, 0,
                dx0/(dx1*(dx0 + dx1)), 0, 0, 0, 0, 0, 0},
            {-((dx1*dy1)/(dx0*(dx0 + dx1)*dy0*(dy0 + dy1))),
                (dx1*dy0 - dx1*dy1)/(dx0*dx0*dy0*dy1 + dx0*dx1*dy0*dy1),
                -((dx1*dy0)/(dx0*(dx0 + dx1)*dy1*(dy0 + dy1))), 0,
                ((-dx0 + dx1)*dy1)/(dx0*dx1*dy0*(dy0 + dy1)),
                ((dx0 - dx1)*(dy0 - dy1))/(dx0*dx1*dy0*dy1),
                ((-dx0 + dx1)*dy0)/(dx0*dx1*dy1*(dy0 + dy1)), 0,
                (dx0*dy1)/(dx1*(dx0 + dx1)*dy0*(dy0 + dy1)),
                (dx0*(-dy0 + dy1))/(dx1*(dx0 + dx1)*dy0*dy1),
                (dx0*dy0)/(dx1*(dx0 + dx1)*dy1*(dy0 + dy1)), 0, 0, 0, 0, 0},
            {(2*dx1)/(dx0*(dx0 + dx1)*dy0*(dy0 + dy1)),
                (dx1*(2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*(dx0 + dx1)*dy0*dy1*dy1*(dy1 + dy2)),
                -((dx1*(dy0 + dy1 + 2*dy2))/(dx0*(dx0 + dx1)*dy1*(dy0 + dy1)*dy2)),
                dx1/(dx0*(dx0 + dx1)*dy2*(dy1 + dy2)),
                (2*dx0 - 2*dx1)/(dx0*dx1*dy0*dy0 + dx0*dx1*dy0*dy1),
                ((dx0 - dx1)*(2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*dx1*dy0*dy1*dy1*(dy1 + dy2)),
                -(((dx0 - dx1)*(dy0 + dy1 + 2*dy2))/(dx0*dx1*dy1*(dy0 + dy1)*dy2)),
                (dx0 - dx1)/(dx0*dx1*dy1*dy2 + dx0*dx1*dy2*dy2),
                (-2*dx0)/(dx1*(dx0 + dx1)*dy0*(dy0 + dy1)),
                -((dx0*(2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx1*(dx0 + dx1)*dy0*dy1*dy1*(dy1 + dy2))),
                (dx0*(dy0 + dy1 + 2*dy2))/(dx1*(dx0 + dx1)*dy1*(dy0 + dy1)*dy2),
                -(dx0/(dx1*(dx0 + dx1)*dy2*(dy1 + dy2))), 0, 0, 0, 0},
            {-(dx1/(dx0*(dx0 + dx1)*dy0*dy1*(dy0 + dy1))),
                -((dx1*(dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*(dx0 + dx1)*dy0*dy1*dy1*dy1*(dy1 + dy2))),
                (dx1*(dy0 + dy1 + dy2))/(dx0*(dx0 + dx1)*dy1*dy1*(dy0 + dy1)*dy2),
                -(dx1/(dx0*(dx0 + dx1)*dy1*dy2*(dy1 + dy2))),
                -((dx0 - dx1)/(dx0*dx1*dy0*dy0*dy1 + dx0*dx1*dy0*dy1*dy1)),
                -(((dx0 - dx1)*(dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*dx1*dy0*dy1*dy1*dy1*(dy1 + dy2))),
                ((dx0 - dx1)*(dy0 + dy1 + dy2))/(dx0*dx1*dy1*dy1*(dy0 + dy1)*dy2),
                -((dx0 - dx1)/(dx0*dx1*dy1*dy1*dy2 + dx0*dx1*dy1*dy2*dy2)),
                dx0/(dx1*(dx0 + dx1)*dy0*dy1*(dy0 + dy1)),
                (dx0*(dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx1*(dx0 + dx1)*dy0*dy1*dy1*dy1*(dy1 + dy2)),
                -((dx0*(dy0 + dy1 + dy2))/(dx1*(dx0 + dx1)*dy1*dy1*(dy0 + dy1)*dy2)),
                dx0/(dx1*(dx0 + dx1)*dy1*dy2*(dy1 + dy2)), 0, 0, 0, 0},
            {0, 2/(dx0*(dx0 + dx1)), 0, 0, 0, -((2/dx0 + 1/(dx1 + dx2))/dx1), 0,
                0, 0, (dx0 + dx1 + 2*dx2)/(dx0*dx1*dx2 + dx1*dx1*dx2), 0, 0, 0,
                -(1/(dx2*(dx1 + dx2))), 0, 0},
            {(2*dy1)/(dx0*(dx0 + dx1)*dy0*(dy0 + dy1)),
                -((2*dy0 - 2*dy1)/(dx0*dx0*dy0*dy1 + dx0*dx1*dy0*dy1)),
                (2*dy0)/(dx0*(dx0 + dx1)*dy1*(dy0 + dy1)), 0,
                ((-2*dx1*(dx1 + dx2) + dx0*(5*dx1 + 6*dx2))*dy1)/(dx0*dx1*dx1*(dx1 + dx2)*dy0*(dy0 + dy1)),
                ((dx0 + 2*(dx1 + dx2))*(dy0 - dy1))/(dx0*dx1*(dx1 + dx2)*dy0*dy1),
                -(((dx0 + 2*(dx1 + dx2))*dy0)/(dx0*dx1*(dx1 + dx2)*dy1*(dy0 + dy1))), 0,
                ((dx0*(dx1 - 6*dx2) + dx1*(dx1 - 4*dx2))*dy1)/(dx1*dx1*(dx0 + dx1)*dx2*dy0*(dy0 + dy1)),
                -(((dx0 + dx1 + 2*dx2)*(dy0 - dy1))/(dx1*(dx0 + dx1)*dx2*dy0*dy1)),
                ((dx0 + dx1 + 2*dx2)*dy0)/(dx1*(dx0 + dx1)*dx2*dy1*(dy0 + dy1)), 0,
                -(dy1/(dx2*(dx1 + dx2)*dy0*(dy0 + dy1))),
                (dy0 - dy1)/(dx1*dx2*dy0*dy1 + dx2*dx2*dy0*dy1), -(dy0/(dx2*(dx1 + dx2)*dy1*(dy0 + dy1))), 0},
            {-4/(dx0*(dx0 + dx1)*dy0*(dy0 + dy1)),
                (-2*(2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*(dx0 + dx1)*dy0*dy1*dy1*(dy1 + dy2)),
                (2*(dy0 + dy1 + 2*dy2))/(dx0*(dx0 + dx1)*dy1*(dy0 + dy1)*dy2),
                -2/(dx0*(dx0 + dx1)*dy2*(dy1 + dy2)),
                (4*dx1*(dx1 + dx2) - 2*dx0*(5*dx1 + 6*dx2))/(dx0*dx1*dx1*(dx1 + dx2)*dy0*(dy0 + dy1)),
                (dx0*(dx1*dy0*(dy1 - 4*dy2) - 6*dx2*dy0*dy2 + 2*dx1*dy1*(dy1 + dy2))
                 + 2*dx1*(dx1 + dx2)*(2*dy1*(dy1 + dy2)
                                      + dy0*(dy1 + 2*dy2)))/(dx0*dx1*dx1*(dx1 + dx2)*dy0*dy1*dy1*(dy1 + dy2)),
                -(((dx0 + 2*(dx1 + dx2))*(dy0 + dy1 + 2*dy2))/(dx0*dx1*(dx1 + dx2)*dy1*(dy0 + dy1)*dy2)),
                (dx0 + 2*(dx1 + dx2))/(dx0*dx1*(dx1 + dx2)*dy2*(dy1 + dy2)),
                (-2*(dx0*(dx1 - 6*dx2) + dx1*(dx1 -4*dx2)))/(dx1*dx1*(dx0 + dx1)*dx2*dy0*(dy0 + dy1)),
                -((dx0*(-6*dx2*dy0*dy2 + 2*dx1*dy1*(dy1 + dy2) + dx1*dy0*(dy1 + 2*dy2))
                   + dx1*(2*dx2*dy0*(dy1 - dy2)+ 2*dx1*dy1*(dy1 + dy2) + 4*dx2*dy1*(dy1 + dy2) + dx1*dy0*(dy1 + 2*dy2)))/(dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1*dy1*(dy1 + dy2))),
                ((dx0 + dx1 + 2*dx2)*(dy0 + dy1 + 2*dy2))/(dx1*(dx0 + dx1)*dx2*dy1*(dy0 + dy1)*dy2),
                -((dx0 + dx1 + 2*dx2)/(dx1*(dx0 + dx1)*dx2*dy2*(dy1 + dy2))), 2/(dx2*(dx1 + dx2)*dy0*(dy0 + dy1)),
                (2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2))/(dx2*(dx1 + dx2)*dy0*dy1*dy1*(dy1 + dy2)),
                -((dy0 + dy1 + 2*dy2)/(dx2*(dx1 + dx2)*dy1*(dy0 + dy1)*dy2)),
                1/(dx2*(dx1 + dx2)*dy2*(dy1 + dy2))},
            {2/(dx0*(dx0 + dx1)*dy0*dy1*(dy0 + dy1)),
                (2*(dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*(dx0 + dx1)*dy0*dy1*dy1*dy1*(dy1 + dy2)),
                (-2*(dy0 + dy1 + dy2))/(dx0*(dx0 + dx1)*dy1*dy1*(dy0 + dy1)*dy2),
                2/(dx0*(dx0 + dx1)*dy1*dy2*(dy1 + dy2)),
                (-2*dx1*(dx1 + dx2) + dx0*(5*dx1 + 6*dx2))/(dx0*dx1*dx1*(dx1 + dx2)*dy0*dy1*(dy0 + dy1)),
                -((dx0*(dx1*dy0*(dy1 - 4*dy2) - 6*dx2*dy0*dy2 + dx1*dy1*(dy1 + dy2))
                   + 2*dx1*(dx1 + dx2)*(dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2)))/(dx0*dx1*dx1*(dx1 + dx2)*dy0*dy1*dy1*dy1*(dy1 + dy2))),
                ((dx0 + 2*(dx1 + dx2))*(dy0 + dy1 + dy2))/(dx0*dx1*(dx1 + dx2)*dy1*dy1*(dy0 + dy1)*dy2),
                -((dx0 + 2*(dx1 + dx2))/(dx0*dx1*(dx1 + dx2)*dy1*dy2*(dy1 + dy2))),
                (dx0*(dx1 - 6*dx2) + dx1*(dx1 - 4*dx2))/(dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1*(dy0 + dy1)),
                (dx0*(-6*dx2*dy0*dy2 + dx1*dy1*(dy1 + dy2) + dx1*dy0*(dy1 + 2*dy2))
                 + dx1*(2*dx2*dy0*(dy1 - dy2) + dx1*dy1*(dy1 + dy2) + 2*dx2*dy1*(dy1 + dy2) + dx1*dy0*(dy1 + 2*dy2)))/(dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1*dy1*dy1*(dy1 + dy2)),
                -(((dx0 + dx1 + 2*dx2)*(dy0 + dy1 + dy2))/(dx1*(dx0 + dx1)*dx2*dy1*dy1*(dy0 + dy1)*dy2)),
                (dx0 + dx1 + 2*dx2)/(dx1*(dx0 + dx1)*dx2*dy1*dy2*(dy1 + dy2)),
                -(1/(dx2*(dx1 + dx2)*dy0*dy1*(dy0 + dy1))),
                -((dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2))/(dx2*(dx1 + dx2)*dy0*dy1*dy1*dy1*(dy1 + dy2))),
                (dy0 + dy1 + dy2)/(dx2*(dx1 + dx2)*dy1*dy1*(dy0 + dy1)*dy2),
                -(1/(dx2*(dx1 + dx2)*dy1*dy2*(dy1 + dy2)))},
            {0, -(1/(dx0*dx1*(dx0 + dx1))), 0, 0, 0, (dx0 + dx1 + dx2)/(dx0*dx1*dx1*(dx1 + dx2)),
                0, 0, 0, -((dx0 + dx1 + dx2)/(dx1*dx1*(dx0 + dx1)*dx2)), 0, 0, 0,
                1/(dx1*dx1*dx2 + dx1*dx2*dx2), 0, 0},
            {-(dy1/(dx0*dx1*(dx0 + dx1)*dy0*(dy0 + dy1))),
                (dy0 - dy1)/(dx0*dx0*dx1*dy0*dy1 + dx0*dx1*dx1*dy0*dy1),
                -(dy0/(dx0*dx1*(dx0 + dx1)*dy1*(dy0 + dy1))), 0,
                ((dx1*(dx1 + dx2) - dx0*(3*dx1 + 4*dx2))*dy1)/(dx0*dx1*dx1*dx1*(dx1 + dx2)*dy0*(dy0 + dy1)),
                -(((dx0 + dx1 + dx2)*(dy0 - dy1))/(dx0*dx1*dx1*(dx1 + dx2)*dy0*dy1)),
                ((dx0 + dx1 + dx2)*dy0)/(dx0*dx1*dx1*(dx1 + dx2)*dy1*(dy0 + dy1)), 0,
                -(((dx0*(dx1 - 4*dx2) + dx1*(dx1 - 3*dx2))*dy1)/(dx1*dx1*dx1*(dx0 + dx1)*dx2*dy0*(dy0 + dy1))),
                ((dx0 + dx1 + dx2)*(dy0 - dy1))/(dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1),
                -(((dx0 + dx1 + dx2)*dy0)/(dx1*dx1*(dx0 + dx1)*dx2*dy1*(dy0 + dy1))), 0,
                dy1/(dx1*dx2*(dx1 + dx2)*dy0*(dy0 + dy1)),
                -((dy0 - dy1)/(dx1*dx1*dx2*dy0*dy1 + dx1*dx2*dx2*dy0*dy1)),
                dy0/(dx1*dx2*(dx1 + dx2)*dy1*(dy0 + dy1)), 0},
            {2/(dx0*dx1*(dx0 + dx1)*dy0*(dy0 + dy1)),
                (2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2))/(dx0*dx1*(dx0 + dx1)*dy0*dy1*dy1*(dy1 + dy2)),
                -((dy0 + dy1 + 2*dy2)/(dx0*dx1*(dx0 + dx1)*dy1*(dy0 + dy1)*dy2)),
                1/(dx0*dx1*(dx0 + dx1)*dy2*(dy1 + dy2)),
                (-2*dx1*(dx1 + dx2) + dx0*(6*dx1 + 8*dx2))/(dx0*dx1*dx1*dx1*(dx1 + dx2)*dy0*(dy0 + dy1)),
                -((dx0*(dx1*dy0*(dy1 - 2*dy2) - 4*dx2*dy0*dy2 + 2*dx1*dy1*(dy1 + dy2))
                   + dx1*(dx1 + dx2)*(2*dy1*(dy1 + dy2)
                                      + dy0*(dy1 + 2*dy2)))/(dx0*dx1*dx1*dx1*(dx1 + dx2)*dy0*dy1*dy1*(dy1 + dy2))),
                ((dx0 + dx1 + dx2)*(dy0 + dy1 + 2*dy2))/(dx0*dx1*dx1*(dx1 + dx2)*dy1*(dy0 + dy1)*dy2),
                -((dx0 + dx1 + dx2)/(dx0*dx1*dx1*(dx1 + dx2)*dy2*(dy1 + dy2))),
                (2*(dx0*(dx1 - 4*dx2) + dx1*(dx1 - 3*dx2)))/(dx1*dx1*dx1*(dx0 + dx1)*dx2*dy0*(dy0 + dy1)),
                (dx0*(-4*dx2*dy0*dy2 + 2*dx1*dy1*(dy1 + dy2) + dx1*dy0*(dy1 + 2*dy2))
                 + dx1*(dx2*dy0*(dy1 - 2*dy2) + 2*dx1*dy1*(dy1 + dy2) + 2*dx2*dy1*(dy1 + dy2)
                        + dx1*dy0*(dy1 + 2*dy2)))/(dx1*dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1*dy1*(dy1 + dy2)),
                -(((dx0 + dx1 + dx2)*(dy0 + dy1 + 2*dy2))/(dx1*dx1*(dx0 + dx1)*dx2*dy1*(dy0 + dy1)*dy2)),
                (dx0 + dx1 + dx2)/(dx1*dx1*(dx0 + dx1)*dx2*dy2*(dy1 + dy2)),
                -2/(dx1*dx2*(dx1 + dx2)*dy0*(dy0 + dy1)),
                -((2*dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2))/(dx1*dx2*(dx1 + dx2)*dy0*dy1*dy1*(dy1 + dy2))),
                (dy0 + dy1 + 2*dy2)/(dx1*dx2*(dx1 + dx2)*dy1*(dy0 + dy1)*dy2),
                -(1/(dx1*dx2*(dx1 + dx2)*dy2*(dy1 + dy2)))},
            {-(1/(dx0*dx1*(dx0 + dx1)*dy0*dy1*(dy0 + dy1))),
                -((dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2))/(dx0*dx1*(dx0 + dx1)*dy0*dy1*dy1*dy1*(dy1 + dy2))),
                (dy0 + dy1 + dy2)/(dx0*dx1*(dx0 + dx1)*dy1*dy1*(dy0 + dy1)*dy2),
                -(1/(dx0*dx1*(dx0 + dx1)*dy1*dy2*(dy1 + dy2))),
                (dx1*(dx1 + dx2) - dx0*(3*dx1 + 4*dx2))/(dx0*dx1*dx1*dx1*(dx1 + dx2)*dy0*dy1*(dy0 + dy1)),
                (dx0*(dx1*dy0*(dy1 - 2*dy2) - 4*dx2*dy0*dy2 + dx1*dy1*(dy1 + dy2))
                 + dx1*(dx1 + dx2)*(dy1*(dy1 + dy2)
                                    + dy0*(dy1 + 2*dy2)))/(dx0*dx1*dx1*dx1*(dx1 + dx2)*dy0*dy1*dy1*dy1*(dy1 + dy2)),
                -(((dx0 + dx1 + dx2)*(dy0 + dy1 + dy2))/(dx0*dx1*dx1*(dx1 + dx2)*dy1*dy1*(dy0 + dy1)*dy2)),
                (dx0 + dx1 + dx2)/(dx0*dx1*dx1*(dx1 + dx2)*dy1*dy2*(dy1 + dy2)),
                -((dx0*(dx1 - 4*dx2) + dx1*(dx1 - 3*dx2))/(dx1*dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1*(dy0 + dy1))),
                -((dx0*(-4*dx2*dy0*dy2 + dx1*dy1*(dy1 + dy2) + dx1*dy0*(dy1 + 2*dy2))
                   + dx1*(dx2*dy0*(dy1 - 2*dy2) + dx1*dy1*(dy1 + dy2) + dx2*dy1*(dy1 + dy2)
                          + dx1*dy0*(dy1 + 2*dy2)))/(dx1*dx1*dx1*(dx0 + dx1)*dx2*dy0*dy1*dy1*dy1*(dy1 + dy2))),
                ((dx0 + dx1 + dx2)*(dy0 + dy1 + dy2))/(dx1*dx1*(dx0 + dx1)*dx2*dy1*dy1*(dy0 + dy1)*dy2),
                -((dx0 + dx1 + dx2)/(dx1*dx1*(dx0 + dx1)*dx2*dy1*dy2*(dy1 + dy2))),
                1/(dx1*dx2*(dx1 + dx2)*dy0*dy1*(dy0 + dy1)),
                (dy1*(dy1 + dy2) + dy0*(dy1 + 2*dy2))/(dx1*dx2*(dx1 + dx2)*dy0*dy1*dy1*dy1*(dy1 + dy2)),
                -((dy0 + dy1 + dy2)/(dx1*dx2*(dx1 + dx2)*dy1*dy1*(dy0 + dy1)*dy2)),
                1/(dx1*dx2*(dx1 + dx2)*dy1*dy2*(dy1 + dy2))}
        };
        double f[16] = {f00, f01, f02, f03, f10, f11, f12, f13, f20, f21, f22, f23, f30, f31, f32, f33};
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j) {
                a[i][j] = 0.;
                for (int k=0; k<16; ++k)
                    a[i][j] += A[4*i+j][k]*f[k];
            }
		  
        wx[0] = dx0; wx[1] = dx1; wx[2] = dx2;
        wy[0] = dy0; wy[1] = dy1; wy[2] = dy2;
    }
    
    double GeneralBicubic::operator()(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                ans += a[i][j]*mypow(x,i)*mypow(y,j);
        
        return ans;
    }
    
    double GeneralBicubic::Dx(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=1; i<4; ++i)
            for (int j=0; j<4; ++j)
                ans += i*a[i][j]*mypow(x,i-1)*mypow(y,j);
        
        return ans;
    }
    
    double GeneralBicubic::Dy(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=0; i<4; ++i)
            for (int j=1; j<4; ++j)
                ans += j*a[i][j]*mypow(x,i)*mypow(y,j-1);
        
        return ans;
    }
    
    double GeneralBicubic::Dxx(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=2; i<4; ++i)
            for (int j=0; j<4; ++j)
                ans += i*(i-1)*a[i][j]*mypow(x,i-2)*mypow(y,j);
        
        return ans;
    }
    
    double GeneralBicubic::Dxy(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=1; i<4; ++i)
            for (int j=1; j<4; ++j)
                ans += i*j*a[i][j]*mypow(x,i-1)*mypow(y,j-1);
        
        return ans;
    }
    
    double GeneralBicubic::Dyy(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=0; i<4; ++i)
            for (int j=2; j<4; ++j)
                ans += j*(j-1)*a[i][j]*mypow(x,i)*mypow(y,j-2);
        
        return ans;
    }
    
    double GeneralBicubic::Dxxx(const double x, const double y) const
    {
        double ans = 0;
        
        for (int j=0; j<4; ++j)
            ans += 6*a[3][j]*mypow(y,j);
        
        return ans;
    }
    
    double GeneralBicubic::Dxxy(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=2; i<4; ++i)
            for (int j=1; j<4; ++j)
                ans += i*(i-1)*j*a[i][j]*mypow(x,i-2)*mypow(y,j-1);
        
        return ans;
    }
    
    double GeneralBicubic::Dxyy(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=1; i<4; ++i)
            for (int j=2; j<4; ++j)
                ans += i*j*(j-1)*a[i][j]*mypow(x,i-1)*mypow(y,j-2);
        
        return ans;
    }
    
    double GeneralBicubic::Dyyy(const double x, const double y) const
    {
        double ans = 0;
        
        for (int i=0; i<4; ++i)
            ans += 6*a[i][3]*mypow(x,i);
        
        return ans;
    }
    
    GeneralBicubic GeneralBicubic::Dx(void)
    {
        GeneralBicubic d;
        
        for (int j=0; j<4; ++j) {
            d.a[3][j] = 0.;
            for (int i=1; i<4; ++i)
                d.a[i-1][j] = i*a[i][j];
        }
        
        return d;
    }
    
    GeneralBicubic GeneralBicubic::Dy(void)
    {
        GeneralBicubic d;
        
        for (int i=0; i<4; ++i) {
            d.a[i][3] = 0.;
            for (int j=1; j<4; ++j)
                d.a[i][j-1] = j*a[i][j];
        }
        
        return d;
    }
    
    inline double sqr(const double x) {return x*x;}
    
    double GeneralBicubic::LocalDist(const double x, const double y,
                              double& ax, double& ay, char& clean) const
    {
        ax = x;
        ay = y;
        double tol = 1.e-3;
        double p, px, py, tax, tay;
        double delta[4] = {1000.,1000.,1000.,1000.};
        int itcount = 0;
        
        while (sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3]))
               > tol*wx[1]*wy[1] && itcount < 20) {
            
            // Find zero
            
            p = F(ax,ay);
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            delta[0] = -p*px/(px*px+py*py);
            delta[1] = -p*py/(px*px+py*py);
            tax = ax+delta[0];
            tay = ay+delta[1];
            
            // Find min distance
            
            F(tax,tay);
            px = Dx(tax,tay);
            py = Dy(tax,tay);
            double denom = px*px+py*py;
            delta[2] = py*(px*(tay-y)-py*(tax-x))/denom;
            delta[3] = -px*(px*(tay-y)-py*(tax-x))/denom;
            
            ax += delta[0]+delta[2];
            ay += delta[1]+delta[3];
            ++itcount;
        }
        
        clean = itcount < 20;
        //clean = itcount < 20 && ax >= 0. && ax <= wx && ay >= 0. && ay <= wy;
        
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }

//ORIGINAL VERSION    
//    double GeneralBicubic::LocalDistPolar(const double r0, const double r, const double theta,
//                                     double& ar, double& atheta, char& clean) const
//    {
//        ar = r;
//        atheta = theta;
//        double tol = 1.e-3;
//        double p, px, py, pr, ptheta, tar, tatheta;
//        double delta[4] = {1000.,1000.,1000.,1000.};
//        int itcount = 0;
//        double x = (r0+r)*cos(theta);
//        double y = (r0+r)*sin(theta);
//        
//        while (sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3]))
//               > tol*wx[1]*wy[1] && itcount < 20) {
//            
//            // Find zero
//            
//            p = F(ar,atheta);
//            pr = Dx(ar,atheta);
//            ptheta = Dy(ar,atheta);
//            px = pr*cos(atheta)-ptheta*sin(atheta)/(r0+ar);
//            py = pr*sin(atheta)+ptheta*cos(atheta)/(r0+ar);
//            delta[0] = -p*px/(px*px+py*py);
//            delta[1] = -p*py/(px*px+py*py);
//            tar = sqrt((r0+ar)*(r0+ar)+2*(r0+ar)*(delta[0]*cos(atheta)+delta[1]*sin(atheta))+delta[0]*delta[0]+delta[1]*delta[1])-r0;
//            tatheta = atan2((r0+ar)*sin(atheta)+delta[1],(r0+ar)*cos(atheta)+delta[0]);
//            delta[0] = tar-ar;
//            delta[1] = tatheta-atheta;
//            
//            // Find min distance
//            
//            p = F(tar,tatheta);
//            pr = Dx(tar,tatheta);
//            ptheta = Dy(tar,tatheta);
//            px = pr*cos(tatheta)-ptheta*sin(tatheta)/(r0+tar);
//            py = pr*sin(tatheta)+ptheta*cos(tatheta)/(r0+tar);
//            double denom = px*px+py*py;
//            delta[2] = py*(px*((r0+tar)*sin(tatheta)-y)-py*((r0+tar)*cos(tatheta)-x))/denom;
//            delta[3] = -px*(px*((r0+tar)*sin(tatheta)-y)-py*((r0+tar)*cos(tatheta)-x))/denom;
//            double tar2 = sqrt((r0+tar)*(r0+tar)+2*(r0+tar)*(delta[2]*cos(tatheta)+delta[3]*sin(tatheta))+delta[2]*delta[2]+delta[3]*delta[3])-r0;
//            double tatheta2 = atan2((r0+tar)*sin(tatheta)+delta[3],(r0+tar)*cos(tatheta)+delta[2]);
//            delta[2] = tar2-tar;
//            delta[3] = tatheta2-tatheta;
//            ar += delta[0]+delta[2];
//            atheta += delta[1]+delta[3];
//            ++itcount;
//        }
//        
//        clean = itcount < 20;
//        //clean = itcount < 20 && ax >= 0. && ax <= wx && ay >= 0. && ay <= wy;
//        
//        return sqrt(sqr(x-(r0+ar)*cos(atheta))+sqr(y-(r0+ar)*sin(atheta)));
//    }

	//MY VERSION
	 double GeneralBicubic::LocalDistPolar(const double r0, const double r, const double theta,
                                     double& ar, double& atheta, char& clean) const
    {
        ar = r;
        atheta = theta;
        double tol = 1.e-3;
        double p, px, py, pr, ptheta, tar, tatheta;
        double delta[4] = {1000.,1000.,1000.,1000.};
        int itcount = 0;
        double x = (r0+r)*cos(theta);
        double y = (r0+r)*sin(theta);
        
        while (sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3]))
               > tol*wx[1]*wy[1] && itcount < 100) {
            
            // Find zero
            p = F(ar,atheta);
            pr = Dx(ar,atheta);
            ptheta = Dy(ar,atheta)/(r0+ar);
				delta[0] = -p*pr/(pr*pr+ptheta*ptheta);
            delta[1] = -p*ptheta/(pr*pr+ptheta*ptheta);
           	tar = ar+delta[0];
				tatheta = atheta+delta[1];

            // Find min distance
           	delta[2] = (r-ar)-pr*((r-ar)*pr+(theta-atheta)*ptheta)/(pr*pr+ptheta*ptheta);
				delta[3] = (theta-atheta)-ptheta*((r-ar)*pr+(theta-atheta)*ptheta)/(pr*pr+ptheta*ptheta);
				
				ar = tar + delta[2];
				atheta = tatheta + delta[3];
            ++itcount;
        }
        
//		  std:: cout << "itcount = " << itcount << std::endl;
        clean = itcount < 100;
        //clean = itcount < 20 && ax >= 0. && ax <= wx && ay >= 0. && ay <= wy;
        
        return sqrt(sqr(x-(r0+ar)*cos(atheta))+sqr(y-(r0+ar)*sin(atheta)));
    }
    
    
    double GeneralBicubic::XCrossing(const double x, const double y, double& ay, char& clean) const
    {
        ay = y;
        double tol = 1.0e-3*wy[1];
        double delta = 1000.;
        int cnt = 0;
        while (fabs(delta) > tol && cnt < 20) {
            delta = -F(x,ay)/Dy(x,ay);
            ay += delta;
            ++cnt;
        }
        clean = fabs(delta) <= tol;
        return ay;
    }
    
    double GeneralBicubic::YCrossing(const double x, const double y, double& ax, char& clean) const
    {
        ax = x;
        double tol = 1.0e-3*wx[1];
        double delta = 1000.;
        int cnt = 0;
        while (fabs(delta) > tol && cnt < 20) {
            delta = -F(ax,y)/Dx(ax,y);
            ax += delta;
            ++cnt;
        }
        clean = fabs(delta) <= tol;
        return ax;
    }
    
}
