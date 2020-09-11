#include <iomanip>
#include "defs.h"
#include "uniformmesh3d.h"
#include "plotwindow3d.h"
#include "plotwindow2d.h"
#include "boundary3d.h"
#include "um3boundary.h"
#ifdef NETCDF
#include "netcdf.hh"
#include "ncvalues.h"
#endif

namespace levelset {
        
    enum {
        LYLZ=0x0001, HYLZ=0x0002, HYHZ=0x0004, LYHZ=0x0008, LXLZ=0x0010,
        HXLZ=0x0020, HXHZ=0x0040, LXHZ=0x0080, LXLY=0x0100, HXLY=0x0200,
        HXHY=0x0400, LXHY=0x0800
    };

#ifndef IO_GLOBALS_DEFINED
#define IO_GLOBALS_DEFINED
    inline int over_value(const double grid_value, const double levelset) 
    {return ((grid_value > levelset) ? 1 : 0);}

#endif

#ifndef NO_GRAPHICS
    //#define DEBUG_IF3D

    void UniformMesh3D::Plot(PlotWindow3D& Port, const int p, 
                             const double lv, const char shaded,
                             const char showboxes) const
    {
        double rng[3][2];
        double vw[3];

        Port.GetRange(rng[0][0], rng[0][1], rng[1][0], rng[1][1], rng[2][0],
                      rng[2][1]);
        Port.GetView(vw[0], vw[1], vw[2]);

        int first[3], last[3], inc[3];

        for (int i=0; i<3; ++i) {
            if (vw[i] < (rng[i][0] + rng[i][1])/2.) {
                switch(i) {
                case 0: first[i] = maxi-2; break;
                case 1: first[i] = maxj-2; break;
                case 2: first[i] = maxk-2; break;
                }
                last[i] = 0;
                inc[i] = -1;
            }
            else {
                switch(i) {
                case 0: last[i] = maxi-2; break;
                case 1: last[i] = maxj-2; break;
                case 2: last[i] = maxk-2; break;
                }
                first[i] = 0;
                inc[i] = 1;
            }
        }

        int box_value, ii, jj, kk;
        static int lookup[127][5] = {1,  273,    0,    0,    0,
                                     1,  392,    0,    0,    0,
                                     2,  152,   25,    0,    0,
                                     1, 2066,    0,    0,    0,
                                     2, 2305, 2051,    0,    0,
                                     2,  392, 2066,    0,    0,
                                     3,  137, 2177, 2051,    0, 
                                     1, 2180,    0,    0,    0,
                                     2,  273, 2180,    0,    0,
                                     2, 2312, 2060,    0,    0,
                                     3, 2060, 2072,   25,    0,
                                     2,  148,   22,    0,    0,
                                     3,  134,  386,  259,    0,
                                     3,  268,  276,   22,    0,
                                     2,   14,   11,    0,    0,
                                     1,  545,    0,    0,    0,
                                     2,  800,  304,    0,    0,
                                     2,  392,  545,    0,    0,
                                     3,  152,  536,  560,    0,
                                     2, 2066,  545,    0,    0,
                                     3,  800,  290, 2306,    0,
                                     3,  392, 2066,  545,    0,
                                     4, 2184, 2568, 2562,  546,
                                     2, 2180,  545,    0,    0,
                                     3, 2180,  800,  304,    0,
                                     3, 2312, 2060,  545,    0,
                                     4, 2060, 2072,  536,  560,
                                     3,   22,  148,  545,    0,
                                     4,  546,  770,  386,  134,
                                     4,  268,  276,   22,  545,
                                     3,  524,  548,   38,    0,
                                     1,  584,    0,    0,    0,
                                     2,  584,  273,    0,    0,
                                     2,  704,  896,    0,    0,
                                     3,  704,  641,  145,    0,
                                     2, 2066,  584,    0,    0,
                                     3, 2305, 2051,  584,    0,
                                     3,  704,  896, 2066,    0,
                                     4, 2051, 2177,  193,  577,
                                     2, 2180,  584,    0,    0,
                                     3,  273, 2180,  584,    0,
                                     3, 2308,  324,  832,    0,
                                     4, 2065, 2561, 2564,  580,
                                     3,  148,   22,  584,    0,
                                     4,  134,  386,  259,  584,
                                     4,   22,  276,  324,  832,
                                     3,   70,  578,  515,    0,
                                     2,  104,   41,    0,    0,
                                     3,  104,  296,  304,    0,
                                     3,  448,  321,   97,    0,
                                     2,  208,  112,    0,    0,
                                     3,  104,   41, 2066,    0,
                                     4, 2306,  290,  352,  328,
                                     4, 2066,  448,  321,   97,
                                     3, 2240, 2114,   98,    0,
                                     3, 2180,  104,   41,    0,
                                     4, 2180,  104,  296,  304,
                                     4,   97,  321, 2368, 2116,
                                     3, 2096, 2084,  100,    0,
                                     4,  148,   22,  104,   41,
                                     3,  392,   98,   70,    0,
                                     3,  273,   98,   70,    0,
                                     2,   98,   70,    0,    0,
                                     1, 1058,    0,    0,    0,
                                     2,  273, 1058,    0,    0,
                                     2,  392, 1058,    0,    0,
                                     3,  152,   25, 1058,    0,
                                     2, 2096, 3104,    0,    0,
                                     3, 2305, 2081, 3104,    0,
                                     3,  392, 2096, 3104,    0,
                                     4, 3104, 2081, 2057, 2184,
                                     2, 2180, 1058,    0,    0,
                                     3,  273, 2180, 1058,    0,
                                     3, 2312, 2060, 1058,    0,
                                     4, 2060, 2072,   25, 1058,
                                     3, 1072, 1044,  148,    0,
                                     4,  289,  416,  164, 1060,
                                     4,  268,  276,   52, 1060,
                                     3, 1036, 1064,   41,    0,
                                     2,  515, 1538,    0,    0,
                                     3,  784,  530, 1538,    0,
                                     3,  392,  515, 1538,    0,
                                     4, 1538,  530,  656,  648,
                                     3, 1537, 1041, 3088,    0,
                                     2, 2816, 3584,    0,    0,
                                     4,  392, 1537, 1041, 3088,
                                     3, 1544, 1160, 3200,    0,
                                     3,  515, 1538, 2180,    0,
                                     4, 2180,  784,  530, 1538,
                                     4,  515, 1538, 2312, 2060,
                                     3, 2066,  524, 1540,    0,
                                     4, 1537, 1041, 1168, 1156,
                                     3,  896,  644, 1540,    0,
                                     3,  273,  524, 1540,    0,
                                     2,  524, 1540,    0,    0,
                                     2,  584, 1058,    0,    0,
                                     3,  584, 1058,  273,    0,
                                     3,  704,  896, 1058,    0,
                                     4, 1058,  704,  641,  145,
                                     3,  584, 2096, 3104,    0,
                                     4,  584, 2305, 2081, 3104,
                                     4,  704,  896, 2096, 3104,
                                     3,  545, 2240, 3136,    0,
                                     3, 2180,  584, 1058,    0,
                                     4,  273, 2180,  584, 1058,
                                     4, 2308,  324,  832, 1058,
                                     3, 2066,  545, 1092,    0,
                                     4,  584, 1072, 1044,  148,
                                     3,  392,  545, 1092,    0,
                                     3,  800,  304, 1092,    0,
                                     2,  545, 1092,    0,    0,
                                     3,   73, 1089, 1027,    0,
                                     4, 1096, 1288, 1282,  274,
                                     4, 1027, 1089,  193,  385,
                                     3, 1216, 1154,  146,    0,
                                     4, 3088, 1041, 1033, 1096,
                                     3, 2312, 2120, 3136,    0,
                                     3,  273, 2240, 3136,    0,
                                     2, 2240, 3136,    0,    0,
                                     4,  137, 2177, 2051, 1092,
                                     3,  392, 2066, 1092,    0,
                                     3, 2305, 2051, 1092,    0,
                                     2, 2066, 1092,    0,    0,
                                     3,  152,   25, 1092,    0,
                                     2,  392, 1092,    0,    0,
                                     2,  273, 1092,    0,    0,
                                     1, 1092,    0,    0,    0};
   
        // assume the port is set up, cleared, and ready to plot

#ifdef DEBUG_IF3D
        ofstream delme("deleteme.out");
#endif
        int maxl = maxi+maxj+maxk-6;
        for (int l=0; l<maxl; ++l) {
            for (int i=0; (i<=l) && (i<maxi-1); ++i) {
                ii = first[0] + i*inc[0];
                for (int j=max(0,l-(2*maxk-2)-i); (j<=l-i) && (j<maxj-1); ++j) {
                    jj = first[1] + j*inc[1];
                    if (l-i-j < maxk-1) {
                        kk = first[2] + (l-i-j)*inc[2];
               
                        if (showboxes) {
                            Port.Line(X(ii), Y(jj), Z(kk), X(ii+1), Y(jj), Z(kk));
                            Port.Line(X(ii+1), Y(jj), Z(kk), X(ii+1), Y(jj+1), Z(kk));
                            Port.Line(X(ii), Y(jj+1), Z(kk), X(ii+1), Y(jj+1), Z(kk));
                            Port.Line(X(ii), Y(jj), Z(kk), X(ii), Y(jj+1), Z(kk));
                            Port.Line(X(ii), Y(jj), Z(kk+1), X(ii+1), Y(jj), Z(kk+1));
                            Port.Line(X(ii+1), Y(jj), Z(kk+1),
                                      X(ii+1), Y(jj+1), Z(kk+1));
                            Port.Line(X(ii), Y(jj+1), Z(kk+1),
                                      X(ii+1), Y(jj+1), Z(kk+1));
                            Port.Line(X(ii), Y(jj), Z(kk+1), X(ii), Y(jj+1), Z(kk+1));
                            Port.Line(X(ii), Y(jj), Z(kk), X(ii), Y(jj), Z(kk+1));
                            Port.Line(X(ii+1), Y(jj), Z(kk), X(ii+1), Y(jj), Z(kk+1));
                            Port.Line(X(ii), Y(jj+1), Z(kk), X(ii), Y(jj+1), Z(kk+1));
                            Port.Line(X(ii+1), Y(jj+1), Z(kk),
                                      X(ii+1), Y(jj+1), Z(kk+1));
                        };

                        box_value = over_value(data_(ii+1,jj+1,kk+1,p), lv);
                        box_value = 2*box_value +over_value(data_(ii+1,jj+1,kk,p), lv);
                        box_value = 2*box_value +over_value(data_(ii+1,jj,kk+1,p),lv);
                        box_value = 2*box_value +over_value(data_(ii+1,jj,kk,p), lv);
                        box_value = 2*box_value +over_value(data_(ii,jj+1,kk+1,p), lv);
                        box_value = 2*box_value +over_value(data_(ii,jj+1,kk,p), lv);
                        box_value = 2*box_value +over_value(data_(ii,jj,kk+1,p),lv);
                        box_value = 2*box_value +over_value(data_(ii,jj,kk,p), lv);
                        if (box_value > 127) box_value = 255 - box_value;

                        double pt[3][3];
               
                        if (box_value != 0) {
                            for (int m=1; m<=lookup[box_value-1][0]; ++m) {
                                int n;
                                for (n=0; n<3; ++n) {
                                    pt[n][0] = X(ii);
                                    pt[n][1] = Y(jj);
                                    pt[n][2] = Z(kk);
                                }
                                n = 0;
                                if (lookup[box_value-1][m] & LYLZ) 
                                    pt[n++][0] += dx*frac(lv, data_(ii,jj,kk,p),
                                                          data_(ii+1,jj,kk,p));
                                if (lookup[box_value-1][m] & HYLZ) {
                                    pt[n][1] += dy;
                                    pt[n++][0] += dx*frac(lv, data_(ii,jj+1,kk,p),
                                                          data_(ii+1,jj+1,kk,p));
                                };
                                if (lookup[box_value-1][m] & HYHZ) {
                                    pt[n][1] += dy;
                                    pt[n][2] += dz;
                                    pt[n++][0] += dx*frac(lv, data_(ii,jj+1,kk+1,p),
                                                          data_(ii+1,jj+1,kk+1,p));
                                };
                                if (lookup[box_value-1][m] & LYHZ) {
                                    pt[n][2] += dz;
                                    pt[n++][0] += dx*frac(lv, data_(ii,jj,kk+1,p),
                                                          data_(ii+1,jj,kk+1,p));
                                };
                                if (lookup[box_value-1][m] & LXLZ) 
                                    pt[n++][1] += dy*frac(lv, data_(ii,jj,kk,p),
                                                          data_(ii,jj+1,kk,p));
                                if (lookup[box_value-1][m] & HXLZ) {
                                    pt[n][0] += dx;
                                    pt[n++][1] += dy*frac(lv, data_(ii+1,jj,kk,p),
                                                          data_(ii+1,jj+1,kk,p));
                                };
                                if (lookup[box_value-1][m] & HXHZ) {
                                    pt[n][0] += dx;
                                    pt[n][2] += dz;
                                    pt[n++][1] += dy*frac(lv, data_(ii+1,jj,kk+1,p),
                                                          data_(ii+1,jj+1,kk+1,p));
                                };
                                if (lookup[box_value-1][m] & LXHZ) {
                                    pt[n][2] += dz;
                                    pt[n++][1] += dy*frac(lv, data_(ii,jj,kk+1,p),
                                                          data_(ii,jj+1,kk+1,p));
                                };
                                if (lookup[box_value-1][m] & LXLY) 
                                    pt[n++][2] += dz*frac(lv, data_(ii,jj,kk,p),
                                                          data_(ii,jj,kk+1,p));
                                if (lookup[box_value-1][m] & HXLY) {
                                    pt[n][0] += dx;
                                    pt[n++][2] += dz*frac(lv, data_(ii+1,jj,kk,p),
                                                          data_(ii+1,jj,kk+1,p));
                                };
                                if (lookup[box_value-1][m] & HXHY) {
                                    pt[n][0] += dx;
                                    pt[n][1] += dy;
                                    pt[n++][2] += dz*frac(lv, data_(ii+1,jj+1,kk,p),
                                                          data_(ii+1,jj+1,kk+1,p));
                                };
                                if (lookup[box_value-1][m] & LXHY) {
                                    pt[n][1] += dy;
                                    pt[n++][2] += dz*frac(lv, data_(ii,jj+1,kk,p),
                                                          data_(ii,jj+1,kk+1,p));
                                };
                                if (n == 3) {
#ifdef DEBUG_IF3D
                                    delme << '(' << ii << ',' << jj << ',' << kk 
                                          << "): (" << pt[0][0] << ',' << pt[0][1] << ','
                                          << pt[0][2] << "), (" << pt[1][0] << ',' 
                                          << pt[1][1] << ',' << pt[1][2] << "), ("
                                          << pt[2][0] << ',' << pt[2][1] << ',' 
                                          << pt[2][2] << ")\n";
#endif
                                    if (shaded) {
                                        Port.Triangle(pt[0][0], pt[0][1], pt[0][2],
                                                      pt[1][0], pt[1][1], pt[1][2],
                                                      pt[2][0], pt[2][1], pt[2][2]);
                                    } else {
                                        Port.Line(pt[0][0], pt[0][1], pt[0][2],
                                                  pt[1][0], pt[1][1], pt[1][2]);
                                        Port.Line(pt[1][0], pt[1][1], pt[1][2],
                                                  pt[2][0], pt[2][1], pt[2][2]);
                                        Port.Line(pt[0][0], pt[0][1], pt[0][2],
                                                  pt[2][0], pt[2][1], pt[2][2]);
                                    }
                                } else {
                                    std::cout << "ERROR: faulty data table for "
                                              << box_value << '\n';
                                };
                            }
                        }
                    }
                }
            }
        }
    }

    void UniformMesh3D::ArrowPlot(PlotWindow3D& Port, const int kx, const int ky,
                                  const int kz, const double s) const
    {
        // assume the port is set up, cleared, and ready to plot
   
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j)
                for (int k=0; k<maxk-1; ++k)
                    Port.Arrow(X(i), Y(j), Z(k), X(i)+s*data_(i,j,k,kx),
                               Y(j)+s*data_(i,j,k,ky), Z(k)+s*data_(i,j,k,kz));
    }
   
    void UniformMesh3D::ArrowPlot(PlotWindow3D& Port, const int kx, const int ky,
                                  const int kz, const int km, const int kd,
                                  const double s) const
    {
        // assume the port is set up, cleared, and ready to plot

        if (km != -1) {
            if (kd != -1) {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j)
                        for (int k=0; k<maxk-1; ++k)
                            Port.Arrow(X(i),Y(j),Z(k),
                                       X(i)+s*data_(i,j,k,kx)*data_(i,j,k,km)
                                       /data_(i,j,k,kd),
                                       Y(j)+s*data_(i,j,k,ky)*data_(i,j,k,km)
                                       /data_(i,j,k,kd),
                                       Z(k)+s*data_(i,j,k,kz)*data_(i,j,k,km)
                                       /data_(i,j,k,kd));
            }
            else {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j)
                        for (int k=0; k<maxk-1; ++k)
                            Port.Arrow(X(i),Y(j),Z(k),
                                       X(i)+s*data_(i,j,k,kx)*data_(i,j,k,km),
                                       Y(j)+s*data_(i,j,k,ky)*data_(i,j,k,km),
                                       Z(k)+s*data_(i,j,k,kz)*data_(i,j,k,km));
            }
        }
        else {
            if (kd != -1) {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j)
                        for (int k=0; k<maxk-1; ++k)
                            Port.Arrow(X(i),Y(j),Z(k),
                                       X(i)+s*data_(i,j,k,kx)/data_(i,j,k,kd),
                                       Y(j)+s*data_(i,j,k,ky)/data_(i,j,k,kd),
                                       Z(k)+s*data_(i,j,k,kz)/data_(i,j,k,kd));
            }
            else {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j)
                        for (int k=0; k<maxk-1; ++k)
                            Port.Arrow(X(i),Y(j),Z(k),
                                       X(i)+s*data_(i,j,k,kx),
                                       Y(j)+s*data_(i,j,k,ky),
                                       Z(k)+s*data_(i,j,k,kz));
            }
        }
         
    }

    void UniformMesh3D::GraphXSlice(PlotWindow3D& Port, const int i, const int l,
                                    const double mult, const char shaded) const
    {
        if (shaded) {
            int j0, jN, jp, k0, kN, kp;
            double view[3];
            Port.GetView(view[0], view[1], view[2]);
            if (view[0] < 0) {
                j0 = maxj-2;
                jN = -1;
                jp = -1;
            } else {
                j0 = 0;
                jN = maxj-1;
                jp = 1;
            }
            if (view[1] < 0) {
                k0 = maxk-2;
                kN = -1;
                kp = -1;
            } else {
                k0 = 0;
                kN = maxk-1;
                kp = 1;
            }
            for (int j=j0; j!=jN; j+=jp)
                for (int k=k0; k!=kN; k+=kp) {
                    Port.Triangle(Y(j),Z(k),data_(i,j,k,l)*mult,Y(j+1),Z(k),data_(i,j+1,k,l)*mult,
                                  Y(j),Z(k+1),data_(i,j,k+1,l)*mult);
                    Port.Triangle(Y(j+1),Z(k),data_(i,j+1,k,l)*mult,Y(j+1),Z(k+1),data_(i,j+1,k+1,l)*mult,
                                  Y(j),Z(k+1),data_(i,j,k+1,l)*mult);
                };
        } else {
            for (int j=0; j<maxj-1; ++j)
                for (int k=0; k<maxk-1; ++k) {
                    Port.Line(Y(j),Z(k),data_(i,j,k,l)*mult,Y(j+1),Z(k),data_(i,j+1,k,l)*mult);
                    Port.Line(Y(j+1),Z(k),data_(i,j+1,k,l)*mult,Y(j+1),Z(k+1),data_(i,j+1,k+1,l)*mult);
                    Port.Line(Y(j),Z(k+1),data_(i,j,k+1,l)*mult,Y(j+1),Z(k+1),data_(i,j+1,k+1,l)*mult);
                    Port.Line(Y(j),Z(k),data_(i,j,k,l)*mult,Y(j),Z(k+1),data_(i,j,k+1,l)*mult);
                };
        }
    }
   
    void UniformMesh3D::GraphYSlice(PlotWindow3D& Port, const int j, const int l,
                                    const double mult, const char shaded) const
    {
        if (shaded) {
            int i0, iN, ip, k0, kN, kp;
            double view[3];
            Port.GetView(view[0], view[1], view[2]);
            if (view[0] < 0) {
                i0 = maxi-2;
                iN = -1;
                ip = -1;
            } else {
                i0 = 0;
                iN = maxi-1;
                ip = 1;
            }
            if (view[1] < 0) {
                k0 = maxk-2;
                kN = -1;
                kp = -1;
            } else {
                k0 = 0;
                kN = maxk-1;
                kp = 1;
            }
            for (int i=i0; i!=iN; i+=ip)
                for (int k=k0; k!=kN; k+=kp) {
                    Port.Triangle(X(i),Z(k),data_(i,j,k,l)*mult,X(i+1),Z(k),data_(i+1,j,k,l)*mult,
                                  X(i),Z(k+1),data_(i,j,k+1,l)*mult);
                    Port.Triangle(X(i+1),Z(k),data_(i+1,j,k,l)*mult,X(i+1),Z(k+1),data_(i+1,j,k+1,l)*mult,
                                  X(i),Z(k+1),data_(i,j,k+1,l)*mult);
                };
        } else {
            for (int i=0; i<maxi-1; ++i)
                for (int k=0; k<maxk-1; ++k) {
                    Port.Line(X(i),Z(k),data_(i,j,k,l)*mult,X(i+1),Z(k),data_(i+1,j,k,l)*mult);
                    Port.Line(X(i+1),Z(k),data_(i+1,j,k,l)*mult,X(i+1),Z(k+1),data_(i+1,j,k+1,l)*mult);
                    Port.Line(X(i),Z(k+1),data_(i,j,k+1,l)*mult,X(i+1),Z(k+1),data_(i+1,j,k+1,l)*mult);
                    Port.Line(X(i),Z(k),data_(i,j,k,l)*mult,X(i),Z(k+1),data_(i,j,k+1,l)*mult);
                };
        }
    }
   
    void UniformMesh3D::GraphZSlice(PlotWindow3D& Port, const int k, const int l,
                                    const double mult, const char shaded) const
    {
        if (shaded) {
            int i0, iN, ip, j0, jN, jp;
            double view[3];
            Port.GetView(view[0], view[1], view[2]);
            if (view[0] < 0) {
                i0 = maxi-2;
                iN = -1;
                ip = -1;
            } else {
                i0 = 0;
                iN = maxi-1;
                ip = 1;
            }
            if (view[1] < 0) {
                j0 = maxj-2;
                jN = -1;
                jp = -1;
            } else {
                j0 = 0;
                jN = maxj-1;
                jp = 1;
            }
            for (int i=i0; i!=iN; i+=ip)
                for (int j=j0; j!=jN; j+=jp) {
                    Port.Triangle(X(i),Y(j),data_(i,j,k,l)*mult,X(i+1),Y(j),data_(i+1,j,k,l)*mult,
                                  X(i),Y(j+1),data_(i,j+1,k,l)*mult);
                    Port.Triangle(X(i+1),Y(j),data_(i+1,j,k,l)*mult,X(i+1),Y(j+1),data_(i+1,j+1,k,l)*mult,
                                  X(i),Y(j+1),data_(i,j+1,k,l)*mult);
                };
        } else {
            for (int i=0; i<maxi-1; ++i)
                for (int j=0; j<maxj-1; ++j) {
                    Port.Line(X(i),Y(j),data_(i,j,k,l)*mult,X(i+1),Y(j),data_(i+1,j,k,l)*mult);
                    Port.Line(X(i+1),Y(j),data_(i+1,j,k,l)*mult,X(i+1),Y(j+1),data_(i+1,j+1,k,l)*mult);
                    Port.Line(X(i),Y(j+1),data_(i,j+1,k,l)*mult,X(i+1),Y(j+1),data_(i+1,j+1,k,l)*mult);
                    Port.Line(X(i),Y(j),data_(i,j,k,l)*mult,X(i),Y(j+1),data_(i,j+1,k,l)*mult);
                };
        }
    }
   
    void UniformMesh3D::ScaledGraphXSlice(PlotWindow3D& Port, const int i, const int l, 
                                          const char shaded) const
    {
        double rng[3][2];
        Port.GetActualRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],rng[2][0],rng[2][1]);
        double mn,mx;
        MinMax(l,mn,mx);
        double mult = (rng[2][1]-rng[2][0])/(mx-mn);
        Port.GetRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],rng[2][0],rng[2][1]);
        Port.SetRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],mult*mn,mult*mx);
        GraphXSlice(Port, i, l, mult, shaded);
        Port.SetRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],rng[2][0],rng[2][1]);
    }
        
    void UniformMesh3D::ScaledGraphYSlice(PlotWindow3D& Port, const int j, const int l, 
                                          const char shaded) const
    {
        double rng[3][2];
        Port.GetActualRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],rng[2][0],rng[2][1]);
        double mn,mx;
        MinMax(l,mn,mx);
        double m1 = rng[2][0]/mn;
        double m2 = rng[2][1]/mx;
        double mult = m1 > 0 ? (m2 > 0 ? min(m1,m2) : m1) : (m2 > 0 ? m2 : 1.);
        GraphYSlice(Port, j, l, mult, shaded);
    }
        
    void UniformMesh3D::ScaledGraphZSlice(PlotWindow3D& Port, const int k, const int l, 
                                          const char shaded) const
    {
        double rng[3][2];
        Port.GetActualRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],rng[2][0],rng[2][1]);
        double mn,mx;
        MinMax(l,mn,mx);
        double m1 = rng[2][0]/mn;
        double m2 = rng[2][1]/mx;
        double mult = m1 > 0 ? (m2 > 0 ? min(m1,m2) : m1) : (m2 > 0 ? m2 : 1.);
        GraphZSlice(Port, k, l, mult, shaded);
    }
        
#endif

    std::ostream& operator<<(std::ostream& s, UniformMesh3D& m)
    {
        s << m.maxi << std::endl << m.maxj << std::endl << m.maxk << std::endl << m.maxl << std::endl;
        s << std::setprecision(20) << m.dx << std::endl
          << std::setprecision(20) << m.dy << std::endl
          << std::setprecision(20) << m.dz << std::endl;
        s << std::setprecision(20) << m.zero[0] << std::endl
          << std::setprecision(20) << m.zero[1] << std::endl
          << std::setprecision(20) << m.zero[2] << std::endl;
        int i;
        for (i=0; i<m.tmaxi*m.tmaxj*m.tmaxk*m.maxl; ++i) s << m.thedata[i];
        return s;
    }
   
    std::istream& operator>>(std::istream& s, UniformMesh3D& m)
    {
        s >> m.tmaxi >> m.tmaxj >> m.maxk >> m.maxl;
        s >> m.dx >> m.dy >> m.dz;
        s >> m.zero[0] >> m.zero[1] >> m.zero[2];
        if (m.thedata) {
            delete[] m.thedata;
            m.BuildWorkGrid(m.maxl);
        }
        for (int i=0; i<m.tmaxi*m.tmaxj*m.tmaxk*m.maxl; ++i) s >> m.thedata[i];
        return s;
    }

#ifdef NETCDF
    void UniformMesh3D::NetCDF(const char* fname, const int kstart,
                               const int kend, const char* comment)
    {
        //cerr << "fname = " << fname << endl;
        int kval[2];
        kval[0] = kstart;
        kval[1] = kend;
        if (kval[0] == -1) {
            kval[0] = 0;
            kval[1] = maxl-1;
        }
        if (kval[1] == -1) kval[1] = kval[0];
        NcFile nfile(fname, NcFile::Replace);
        if (comment != NULL) 
            nfile.add_att("comment", comment);
        NcDim* xdim = nfile.add_dim("x", tmaxi);
        NcDim* ydim = nfile.add_dim("y", tmaxj);
        NcDim* zdim = nfile.add_dim("z", tmaxk);
        NcDim* kdim = nfile.add_dim("k", kval[1]-kval[0]+1);
        NcVar* var = nfile.add_var("data", ncDouble, kdim, xdim, ydim, zdim);
        NcVar* xcoord = nfile.add_var("xcoords", ncDouble, xdim);
        NcVar* ycoord = nfile.add_var("ycoords", ncDouble, ydim);
        NcVar* zcoord = nfile.add_var("zcoords", ncDouble, zdim);
        var->put(&(thedata[kval[0]*tmaxi*tmaxj*tmaxk]), kval[1]-kval[0]+1, tmaxi, tmaxj, tmaxk);
        double* coords = new double[tmaxi+tmaxj+tmaxk];
        int xlo = bc->Width(BDRY3D_XLO);
        int ylo = bc->Width(BDRY3D_YLO);
        int zlo = bc->Width(BDRY3D_ZLO);
        for (int k=0; k<tmaxk; ++k) coords[k] = Z(k-zlo);   
        for (int j=0; j<tmaxj; ++j) coords[tmaxk+j] = Y(j-ylo);
        for (int i=0; i<tmaxi; ++i) coords[tmaxk+tmaxj+i] = X(i-xlo);
        xcoord->put(&(coords[tmaxj+tmaxk]), tmaxi);
        ycoord->put(&(coords[tmaxk]), tmaxj);
        zcoord->put(coords, tmaxk);
        delete[] coords;
    }

    void UniformMesh3D::ReadNetCDF(const char* fname, const int kstart)
    {
        NcFile nfile(fname);
        NcDim* kdim = nfile.get_dim("k"); // number of layers
        NcVar* var = nfile.get_var("data");
        var->get(&(thedata[kstart*tmaxi*tmaxj*tmaxk]), kdim->size(), tmaxi, tmaxj, tmaxk);
    }
#endif

    void UniformMesh3D::WriteBinary(std::string fname, const int kstart,
                                    const int kend, const char* comment) const
    {
        WriteBinary(fname.c_str(),kstart,kend,comment);
    }

    void UniformMesh3D::WriteBinary(const char* fname, const int kstart,
                                    const int kend, const char* comment) const
    {
        int kval[2];
        kval[0] = kstart;
        kval[1] = kend;
        if (kval[0] == -1) {
            kval[0] = 0;
            kval[1] = maxl-1;
        }
        if (kval[1] == -1) kval[1] = kval[0];

        std::ofstream out(fname, std::ios::out | std::ios::binary);
        out.write(reinterpret_cast<const char*>(&tmaxi), sizeof(tmaxi));
        out.write(reinterpret_cast<const char*>(&tmaxi), sizeof(tmaxj));
        out.write(reinterpret_cast<const char*>(&tmaxi), sizeof(tmaxk));
        int klen = kval[1]-kval[0]+1;
        out.write(reinterpret_cast<char*>(&klen),sizeof(klen));
        out.write(reinterpret_cast<char*>(&(thedata[kval[0]*tmaxi*tmaxj*tmaxk])),
                  klen*tmaxi*tmaxj*tmaxk*sizeof(double));
        double* coords = new double[tmaxi+tmaxj+tmaxk];
        int xlo = bc->Width(BDRY3D_XLO);
        int ylo = bc->Width(BDRY3D_YLO);
        int zlo = bc->Width(BDRY3D_ZLO);
        for (int i=0; i<tmaxi; ++i) coords[i] = X(i-xlo);
        for (int j=0; j<tmaxj; ++j) coords[tmaxi+j] = Y(j-ylo);
        for (int k=0; k<tmaxk; ++k) coords[tmaxi+tmaxj+k] = Z(k-zlo);
        out.write(reinterpret_cast<char*>(coords),(tmaxi+tmaxj+tmaxk)*sizeof(double));
        out.close();
        delete[] coords;
    }
}
