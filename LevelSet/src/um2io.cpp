#include <iomanip>
#ifdef NETCDF
#include "netcdf.hh"
#include "ncvalues.h"
#endif
#include "defs.h"
#include "uniformmesh2d.h"
#include "plotwindow2d.h"
#include "um2boundary.h"

namespace levelset {
        
#ifndef IO_GLOBALS_DEFINED
#define IO_GLOBALS_DEFINED
    inline int over_value(const double grid_value, const double levelset) 
    {return ((grid_value > levelset) ? 1 : 0);}

    inline double fraction_(const double grid_low, const double grid_high,
                            const double value) 
    {return (value-grid_low)/(grid_high-grid_low);}
#endif

#ifndef NO_GRAPHICS

    void UniformMesh2D::Plot(PlotWindow2D& Port, const int k, 
                             const double lv, const char showboxes) const
    {
        int box_value;
   
        // assume the port is set up, cleared, and ready to plot
   
        if (showboxes) {
            for (int i=0; i<maxi; ++i)
                Port.Line(X(i), Y(0), X(i), Y(maxj-1));
            for (int j=0; j<maxj; ++j)
                Port.Line(X(0), Y(j), X(maxi-1), Y(j));
        };
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
      
                box_value = over_value(data_(i,j,k), lv);
                box_value = 2*box_value +over_value(data_(i+1,j,k), lv);
                box_value = 2*box_value+over_value(data_(i+1,j+1,k),lv);
                box_value = 2*box_value +over_value(data_(i,j+1,k), lv);
                if (box_value > 7) box_value = 15 - box_value;
      
                switch(box_value) {
                case 1:
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    break;
                case 2:
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    break;
                case 5:
                {
                    double tx[4];
                    double ty[4];
                    double td[3];
                    tx[0] = X(i+fraction_(data_(i,j+1,k),data_(i+1,j+1,k),lv));
                    ty[0] = Y(j+1);
                    tx[1] = X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv));
                    ty[1] = Y(j);
                    tx[2] = X(i);
                    ty[2] = Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv));
                    tx[3] =  X(i+1);
                    ty[3] = Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv));
                    td[0] = sqrt(sqr(tx[0]-tx[1])+sqr(ty[0]-ty[1]))
                        +sqrt(sqr(tx[2]-tx[3])+sqr(ty[2]-ty[3]));
                    td[1] = sqrt(sqr(tx[0]-tx[2])+sqr(ty[0]-ty[2]))
                        +sqrt(sqr(tx[1]-tx[3])+sqr(ty[1]-ty[3]));
                    td[2] = sqrt(sqr(tx[0]-tx[3])+sqr(ty[0]-ty[3]))
                        +sqrt(sqr(tx[2]-tx[1])+sqr(ty[2]-ty[1]));
                    if (td[0] < td[1] && td[0] < td[2]) {
                        Port.Line(tx[0], ty[0], tx[1], ty[1]);
                        Port.Line(tx[2], ty[2], tx[3], ty[3]);
                    }
                    else if (td[1] <= td[0] && td[1] < td[2]) {
                        Port.Line(tx[0], ty[0], tx[2], ty[2]);
                        Port.Line(tx[1], ty[1], tx[3], ty[3]);
                    }
                    else {
                        Port.Line(tx[0], ty[0], tx[3], ty[3]);
                        Port.Line(tx[1], ty[1], tx[2], ty[2]);
                    }
                }
                break;
                case 3:
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)));
                    break;
                case 4:
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 6:
                    Port.Line(X(i+fraction_(data_(i,j+1,k),data_(i+1,j+1,k),lv)), 
                              Y(j+1),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 7:
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)), 
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                };
            };
    }

    void UniformMesh2D::Plot2(PlotWindow2D& Port, const int k, 
                              const double lv, const int reflvl,
                              const char showboxes) const
    {
        int box_value;
        double a[2][2];
        double xx, yy;
        double ldx, ldy;
   
        // assume the port is set up, cleared, and ready to plot
   
        if (showboxes) {
            for (int i=0; i<maxi; ++i)
                Port.Line(X(i), Y(0), X(i), Y(maxj-1));
            for (int j=0; j<maxj; ++j)
                Port.Line(X(0), Y(j), X(maxi-1), Y(j));
        };
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
      
                box_value = over_value(data_(i,j,k), lv);
                box_value = 2*box_value +over_value(data_(i+1,j,k), lv);
                box_value = 2*box_value+over_value(data_(i+1,j+1,k),lv);
                box_value = 2*box_value +over_value(data_(i,j+1,k), lv);
                if (box_value > 7) box_value = 15 - box_value;
                if (box_value > 0) {
                    Bicubic p = Interp2(i,j,k);
                    ldx = dx/reflvl;
                    ldy = dy/reflvl;
                    for (int m=0; m<reflvl; ++m)
                        for (int n=0; n<reflvl; ++n) {
                            xx = X(i)+m*ldx;
                            yy = Y(j)+n*ldy;
                            a[0][0] = p(m*ldx,n*ldy);
                            a[1][0] = p((m+1)*ldx,n*ldy); 
                            a[0][1] = p(m*dx/reflvl,(n+1)*dy/reflvl);
                            a[1][1] = p((m+1)*dx/reflvl,(n+1)*dy/reflvl);
                            box_value = over_value(a[0][0], lv);
                            box_value = 2*box_value +over_value(a[1][0], lv);
                            box_value = 2*box_value+over_value(a[1][1],lv);
                            box_value = 2*box_value +over_value(a[0][1], lv);
                            if (box_value > 7) box_value = 15 - box_value;
                            switch(box_value) {
                            case 1:
                                Port.Line(xx, yy+ldy*fraction_(a[0][0], a[0][1], lv),
                                          xx+ldx*fraction_(a[0][1], a[1][1], lv), yy+ldy);
                                break;
                            case 2:
                                Port.Line(xx+ldx, yy+ldy*fraction_(a[1][0], a[1][1],lv),
                                          xx+ldx*fraction_(a[0][1], a[1][1],lv), yy+ldy);
                                break;
                            case 5:
                            {
                                double tx[4];// DAVE GO HERE
                                double ty[4];
                                double td[3];
                                tx[0] = xx+ldx*fraction_(a[0][1],a[1][1],lv);
                                ty[0] = yy+ldy;
                                tx[1] = xx+ldx*fraction_(a[0][0], a[1][0],lv);
                                ty[1] = yy;
                                tx[2] = xx;
                                ty[2] = yy+ldy*fraction_(a[0][0], a[0][1],lv);
                                tx[3] =  xx+ldx;
                                ty[3] = yy+ldy*fraction_(a[1][0], a[1][1],lv);
                                td[0] = sqrt(sqr(tx[0]-tx[1])+sqr(ty[0]-ty[1]))
                                    +sqrt(sqr(tx[2]-tx[3])+sqr(ty[2]-ty[3]));
                                td[1] = sqrt(sqr(tx[0]-tx[2])+sqr(ty[0]-ty[2]))
                                    +sqrt(sqr(tx[1]-tx[3])+sqr(ty[1]-ty[3]));
                                td[2] = sqrt(sqr(tx[0]-tx[3])+sqr(ty[0]-ty[3]))
                                    +sqrt(sqr(tx[2]-tx[1])+sqr(ty[2]-ty[1]));
                                if (td[0] < td[1] && td[0] < td[2]) {
                                    Port.Line(tx[0], ty[0], tx[1], ty[1]);
                                    Port.Line(tx[2], ty[2], tx[3], ty[3]);
                                }
                                else if (td[1] <= td[0] && td[1] < td[2]) {
                                    Port.Line(tx[0], ty[0], tx[2], ty[2]);
                                    Port.Line(tx[1], ty[1], tx[3], ty[3]);
                                }
                                else {
                                    Port.Line(tx[0], ty[0], tx[3], ty[3]);
                                    Port.Line(tx[1], ty[1], tx[2], ty[2]);
                                }
                            }
                            break;
                            case 3:
                                Port.Line(xx, yy+ldy*fraction_(a[0][0], a[0][1],lv),
                                          xx+ldx, yy+ldy*fraction_(a[1][0], a[1][1],lv));
                                break;
                            case 4:
                                Port.Line(xx+ldx, yy+ldy*fraction_(a[1][0], a[1][1],lv),
                                          xx+ldx*fraction_(a[0][0], a[1][0],lv), yy);
                                break;
                            case 6:
                                Port.Line(xx+ldx*fraction_(a[0][1],a[1][1],lv), yy+ldy,
                                          xx+ldx*fraction_(a[0][0], a[1][0],lv), yy);
                                break;
                            case 7:
                                Port.Line(xx, yy+ldy*fraction_(a[0][0], a[0][1],lv), 
                                          xx+ldx*fraction_(a[0][0], a[1][0],lv), yy);
                                break;
                            };
                        };
                }
            }
    }

    void UniformMesh2D::FilledPlot(PlotWindow2D& Port, const int k, 
                                   const unsigned short shade,
                                   const double lv, const char showboxes) const
    {
        int box_value;
        double xpoly[5], ypoly[5];
   
        // assume the port is set up, cleared, and ready to plot
   
        if (showboxes) {
            for (int i=0; i<maxi; ++i)
                Port.Line(X(i), Y(0), X(i), Y(maxj-1));
            for (int j=0; j<maxj; ++j)
                Port.Line(X(0), Y(j), X(maxi-1), Y(j));
        };
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
      
                box_value = over_value(data_(i,j,k), lv);
                box_value = 2*box_value +over_value(data_(i+1,j,k), lv);
                box_value = 2*box_value+over_value(data_(i+1,j+1,k),lv);
                box_value = 2*box_value +over_value(data_(i,j+1,k), lv);
      
                switch(box_value) {
                case 1:
                    Port.Triangle(X(i), Y(j+1), 
                                  X(i), Y(j+fraction_(data_(i,j,k),data_(i,j+1,k),lv)),
                                  X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                                  Y(j+1), shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    break;
                case 2:
                    Port.Triangle(X(i+1), Y(j+1), X(i+1),
                                  Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                                  X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                                  Y(j+1), shade);
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    break;
                case 3:
                    xpoly[0] = X(i+1); ypoly[0] = Y(j+1);
                    xpoly[1] = X(i); ypoly[1] = Y(j+1);
                    xpoly[2] = X(i); 
                    ypoly[2] = Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv));
                    xpoly[3] = X(i+1);
                    ypoly[3] = Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv));
                    Port.Polygon(4, xpoly, ypoly, shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)));
                    break;
                case 4:
                    Port.Triangle(X(i+1), Y(j), X(i+1), 
                                  Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                                  X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), 
                                  Y(j), shade);
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 5:
                    Port.Triangle(X(i), Y(j+1), 
                                  X(i), Y(j+fraction_(data_(i,j,k),data_(i,j+1,k),lv)),
                                  X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                                  Y(j+1), shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    Port.Triangle(X(i+1), Y(j), X(i+1), 
                                  Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                                  X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), 
                                  Y(j), shade);
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 6:
                    xpoly[0] = X(i+1); ypoly[0] = Y(j);
                    xpoly[1] = X(i+1); ypoly[1] = Y(j+1);
                    xpoly[2] = X(i+fraction_(data_(i,j+1,k),data_(i+1,j+1,k),lv));
                    ypoly[2] = Y(j+1);
                    xpoly[3] = X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv));
                    ypoly[3] = Y(j);
                    Port.Polygon(4, xpoly, ypoly, shade);
                    Port.Line(X(i+fraction_(data_(i,j+1,k),data_(i+1,j+1,k),lv)), 
                              Y(j+1),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 7:
                    xpoly[0] = X(i+1); ypoly[0] = Y(j);
                    xpoly[1] = X(i+1); ypoly[1] = Y(j+1);
                    xpoly[2] = X(i); ypoly[2] = Y(j+1);
                    xpoly[3] = X(i);
                    ypoly[3] = Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv));
                    xpoly[4] = X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv));
                    ypoly[4] = Y(j);
                    Port.Polygon(5, xpoly, ypoly, shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)), 
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 8:
                    Port.Triangle(X(i), Y(j), X(i), 
                                  Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)), 
                                  X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), 
                                  Y(j), shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)), 
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;            
                case 9:
                    xpoly[0] = X(i); ypoly[0] = Y(j);
                    xpoly[1] = X(i); ypoly[1] = Y(j+1);
                    xpoly[2] = X(i+fraction_(data_(i,j+1,k),data_(i+1,j+1,k),lv));
                    ypoly[2] = Y(j+1);
                    xpoly[3] = X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv));
                    ypoly[3] = Y(j);
                    Port.Polygon(4, xpoly, ypoly, shade);
                    Port.Line(X(i+fraction_(data_(i,j+1,k),data_(i+1,j+1,k),lv)), 
                              Y(j+1),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 10:
                    Port.Triangle(X(i+1), Y(j+1), X(i+1),
                                  Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                                  X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                                  Y(j+1), shade);
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    Port.Triangle(X(i), Y(j), X(i), 
                                  Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)), 
                                  X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), 
                                  Y(j), shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)), 
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;            
                case 11:
                    xpoly[0] = X(i); ypoly[0] = Y(j);
                    xpoly[1] = X(i); ypoly[1] = Y(j+1);
                    xpoly[2] = X(i+1); ypoly[2] = Y(j+1);
                    xpoly[3] = X(i+1);
                    ypoly[3] = Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv));
                    xpoly[4] = X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv));
                    ypoly[4] = Y(j);
                    Port.Polygon(5, xpoly, ypoly, shade);
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j,k), data_(i+1,j,k),lv)), Y(j));
                    break;
                case 12:
                    xpoly[0] = X(i+1); ypoly[0] = Y(j);
                    xpoly[1] = X(i); ypoly[1] = Y(j);
                    xpoly[2] = X(i); 
                    ypoly[2] = Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv));
                    xpoly[3] = X(i+1);
                    ypoly[3] = Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv));
                    Port.Polygon(4, xpoly, ypoly, shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)));
                    break;            
                case 13:
                    xpoly[0] = X(i); ypoly[0] = Y(j+1);
                    xpoly[1] = X(i); ypoly[1] = Y(j);
                    xpoly[2] = X(i+1); ypoly[2] = Y(j);
                    xpoly[3] = X(i+1);
                    ypoly[3] = Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv));
                    xpoly[4] = X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv));
                    ypoly[4] = Y(j+1);
                    Port.Polygon(5, xpoly, ypoly, shade);
                    Port.Line(X(i+1), 
                              Y(j+fraction_(data_(i+1,j,k), data_(i+1,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    break;          
                case 14:
                    xpoly[0] = X(i+1); ypoly[0] = Y(j+1);
                    xpoly[1] = X(i+1); ypoly[1] = Y(j);
                    xpoly[2] = X(i); ypoly[2] = Y(j);
                    xpoly[3] = X(i);
                    ypoly[3] = Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv));
                    xpoly[4] = X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv));
                    ypoly[4] = Y(j+1);
                    Port.Polygon(5, xpoly, ypoly, shade);
                    Port.Line(X(i), Y(j+fraction_(data_(i,j,k), data_(i,j+1,k),lv)),
                              X(i+fraction_(data_(i,j+1,k), data_(i+1,j+1,k),lv)),
                              Y(j+1));
                    break;            
                case 15:
                    xpoly[0] = X(i); ypoly[0] = Y(j);
                    xpoly[1] = X(i+1); ypoly[1] = Y(j);
                    xpoly[2] = X(i+1); ypoly[2] = Y(j+1);
                    xpoly[3] = X(i); ypoly[3] = Y(j+1);
                    Port.Polygon(4, xpoly, ypoly, shade);
                    break;
                };
            };
    }

    void UniformMesh2D::ScaledGraph(PlotWindow3D& Port, const int k, const char shaded) const
    {
        double rng[3][2];
        Port.GetActualRange(rng[0][0],rng[0][1],rng[1][0],rng[1][1],rng[2][0],rng[2][1]);
        double mn,mx;
        MinMax(k,mn,mx);
        double m1 = rng[2][0]/mn;
        double m2 = rng[2][1]/mx;
        double mult = m1 > 0 ? (m2 > 0 ? min(m1,m2) : m1) : (m2 > 0 ? m2 : 1.);
        Graph(Port, k, mult, shaded);
    }
        

    void UniformMesh2D::Graph(PlotWindow3D& Port, const int k, const double mult, 
                              const char shaded) const
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
                    Port.Triangle(X(i),Y(j),data_(i,j,k)*mult,X(i+1),Y(j),data_(i+1,j,k)*mult,
                                  X(i),Y(j+1),data_(i,j+1,k)*mult);
                    Port.Triangle(X(i+1),Y(j),data_(i+1,j,k)*mult,X(i+1),Y(j+1),data_(i+1,j+1,k)*mult,
                                  X(i),Y(j+1),data_(i,j+1,k)*mult);
                };
        } else {
            for (int i=0; i<maxi-1; ++i)
                for (int j=0; j<maxj-1; ++j) {
                    Port.Line(X(i),Y(j),data_(i,j,k)*mult,X(i+1),Y(j),data_(i+1,j,k)*mult);
                    Port.Line(X(i+1),Y(j),data_(i+1,j,k)*mult,X(i+1),Y(j+1),data_(i+1,j+1,k)*mult);
                    Port.Line(X(i),Y(j+1),data_(i,j+1,k)*mult,X(i+1),Y(j+1),data_(i+1,j+1,k)*mult);
                    Port.Line(X(i),Y(j),data_(i,j,k)*mult,X(i),Y(j+1),data_(i,j+1,k)*mult);
                };
        }
    }

    void UniformMesh2D::ArrowPlot(PlotWindow2D& Port, const int kx, const int ky,
                                  const double s) const
    {
        // assume the port is set up, cleared, and ready to plot
   
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) 
                Port.Arrow(X(i),Y(j),
                           X(i)+s*data_(i,j,kx),
                           Y(j)+s*data_(i,j,ky));
    }
   
    void UniformMesh2D::ArrowPlot(PlotWindow2D& Port, const int kx, const int ky,
                                  const int km, const int kd, const double s) const
    {
        // assume the port is set up, cleared, and ready to plot

        if (km != -1) {
            if (kd != -1) {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j) 
                        Port.Arrow(X(i),Y(j),
                                   X(i)
                                   +s*data_(i,j,kx)*data_(i,j,km)/data_(i,j,kd),
                                   Y(j)
                                   +s*data_(i,j,ky)*data_(i,j,km)/data_(i,j,kd));
            }
            else {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j) 
                        Port.Arrow(X(i),Y(j),
                                   X(i)+s*data_(i,j,kx)*data_(i,j,km),
                                   Y(j)+s*data_(i,j,ky)*data_(i,j,km));
            }
        }
        else {
            if (kd != -1) {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j) 
                        Port.Arrow(X(i),Y(j),
                                   X(i)+s*data_(i,j,kx)/data_(i,j,kd),
                                   Y(j)+s*data_(i,j,ky)/data_(i,j,kd));
            }
            else {
                for (int i=0; i<maxi-1; ++i)
                    for (int j=0; j<maxj-1; ++j) 
                        Port.Arrow(X(i),Y(j),
                                   X(i)+s*data_(i,j,kx),
                                   Y(j)+s*data_(i,j,ky));
            }
        }
         
    }
   
#endif

    void UniformMesh2D::ShadePlot(PlotWindow2D& Port, const int k, 
                                  const double white, const double black, 
                                  const char showboxes) const
    {
        double x[4], y[4];
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j) {
                double gray = (data_(i,j,k)-black)/(white-black);
                gray = min(gray, 1.);
                gray = max(gray, 0.);
                unsigned short shade = gray*(USHRT_MAX-1);
                x[0] = max(X(0),X(i)-dx/2);
                x[1] = min(X(maxi-1),X(i)+dx/2);
                x[2] = x[1];
                x[3] = x[0];
                y[0] = max(Y(0),Y(j)-dy/2);
                y[1] = y[0];
                y[2] = min(Y(maxj-1),Y(j)+dy/2);
                y[3] = y[2];
                Port.Polygon(4, x, y, shade);
            }
        
        if (showboxes) {
            for (int i=0; i<maxi; ++i)
                Port.Line(X(i), Y(0), X(i), Y(maxj-1));
            for (int j=0; j<maxj; ++j)
                Port.Line(X(0), Y(j), X(maxi-1), Y(j));
        };
    }

#define PORT(I,J) Port[(maxj-1-(J))*(maxi+1)+(I)]

    void UniformMesh2D::Plot(std::ostream& s, const int k, const double lv, 
                             const char showboxes) const
    {
        int         box_value;
        char*    Port;

        // Set up the text port

        Port = new char[(maxi+1)*maxj];
        for (int j=0; j<maxj; ++j) {
            for (int i=0; i<maxi; ++i)
                PORT(i,j) = ' ';
            PORT(maxi,j) = '\n';
        }
        Port[(maxi+1)*maxj-1] = '\0';

        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
                if (showboxes) PORT(i, j) = ':';
   
                box_value = over_value(data_(i,j,k), lv);
                box_value = 2*box_value +over_value(data_(i+1,j,k), lv);
                box_value = 2*box_value+over_value(data_(i+1,j+1,k),lv);
                box_value = 2*box_value + over_value(data_(i,j+1,k),lv);
                if (box_value != 0 && box_value != 15) 
                    PORT(i, j) = '*';
            }
        s << Port << std::endl;
        delete[] Port;
    }

    std::ostream& operator<<(std::ostream& s, UniformMesh2D& m)
    {
        s << m.maxi << std::endl << m.maxj << std::endl << m.maxk << std::endl;
        s << std::setprecision(20) << m.dx << std::endl << std::setprecision(20) << m.dy << std::endl;
        s << std::setprecision(20) << m.zero[0] << std::endl
          << std::setprecision(20) << m.zero[1] << std::endl;
        int i;
        for (i=0; i<m.tmaxi*m.tmaxj*m.maxk; ++i) s << m.thedata[i];
        return s;
    }
   
    std::istream& operator>>(std::istream& s, UniformMesh2D& m)
    {
        s >> m.tmaxi >> m.tmaxj >> m.maxk;
        s >> m.dx >> m.dy;
        s >> m.zero[0] >> m.zero[1];
        if (m.thedata) {
            delete[] m.thedata;
            m.BuildWorkGrid(m.maxk);
        }
        for (int i=0; i<m.tmaxi*m.tmaxj*m.maxk; ++i) s >> m.thedata[i];
        return s;
    }

#ifdef NETCDF
    void UniformMesh2D::NetCDF(const char* fname, const int kstart,
                               const int kend, const char* comment)
    {
        //cerr << "fname = " << fname << endl;
        int kval[2];
        kval[0] = kstart;
        kval[1] = kend;
        if (kval[0] == -1) {
            kval[0] = 0;
            kval[1] = maxk-1;
        }
        if (kval[1] == -1) kval[1] = kval[0];
        NcFile nfile(fname, NcFile::Replace);
        if (comment != NULL) 
            nfile.add_att("comment", comment);
        NcDim* xdim = nfile.add_dim("x", tmaxi);
        NcDim* ydim = nfile.add_dim("y", tmaxj);
        NcDim* kdim = nfile.add_dim("k", kval[1]-kval[0]+1);
        NcVar* var = nfile.add_var("data", ncDouble, kdim, xdim, ydim);
        NcVar* xcoord = nfile.add_var("xcoords", ncDouble, xdim);
        NcVar* ycoord = nfile.add_var("ycoords", ncDouble, ydim);
        var->put(&(thedata[kval[0]*tmaxi*tmaxj]), kval[1]-kval[0]+1, tmaxi, tmaxj);
        double* coords = new double[tmaxi+tmaxj];
        int xlo = bc->Width(BDRY2D_XLO);
        int ylo = bc->Width(BDRY2D_YLO);
        for (int i=0; i<tmaxi; ++i) coords[i] = X(i-xlo);   
        for (int j=0; j<tmaxj; ++j) coords[tmaxi+j] = Y(j-ylo);
        xcoord->put(coords, tmaxi);
        ycoord->put(&(coords[tmaxi]), tmaxj);
        delete[] coords;
    }

    void UniformMesh2D::ReadNetCDF(const char* fname, const int kstart)
    {
        NcFile nfile(fname);
        NcDim* kdim = nfile.get_dim("k"); // number of layers
        NcVar* var = nfile.get_var("data");
        var->get(&(thedata[kstart*tmaxi*tmaxj]), kdim->size(), tmaxi, tmaxj);
    }

#endif

    void UniformMesh2D::WriteBinary(std::string fname, const int kstart,
                                    const int kend, const char* comment)
    {
        WriteBinary(fname.c_str(),kstart,kend,comment);
    }

    void UniformMesh2D::WriteBinary(const char* fname, const int kstart,
                                    const int kend, const char* comment)
    {
        int kval[2];
        kval[0] = kstart;
        kval[1] = kend;
        if (kval[0] == -1) {
            kval[0] = 0;
            kval[1] = maxk-1;
        }
        if (kval[1] == -1) kval[1] = kval[0];

        //std::cout << "Writing: " << fname << std::endl;
   
        std::ofstream out(fname, std::ios::out | std::ios::binary);
        out.write(reinterpret_cast<char*>(&tmaxi),sizeof(tmaxi));
        out.write(reinterpret_cast<char*>(&tmaxj),sizeof(tmaxj));
        int klen = kval[1]-kval[0]+1;
        out.write(reinterpret_cast<char*>(&klen),sizeof(klen));
        out.write(reinterpret_cast<char*>(&(thedata[kval[0]*tmaxi*tmaxj])),
                  klen*tmaxi*tmaxj*sizeof(double));
        double* coords = new double[tmaxi+tmaxj];
        int xlo = bc->Width(BDRY2D_XLO);
        int ylo = bc->Width(BDRY2D_YLO);
        for (int i=0; i<tmaxi; ++i) coords[i] = X(i-xlo);   
        for (int j=0; j<tmaxj; ++j) coords[tmaxi+j] = Y(j-ylo);
        out.write(reinterpret_cast<char*>(coords),(tmaxi+tmaxj)*sizeof(double));
        out.close();
        delete[] coords;
    }

    void UniformMesh2D::ReadBinary(std::string fname, const int kstart)
    {
        ReadBinary(fname.c_str(),kstart);
    }

    void UniformMesh2D::ReadBinary(const char* fname, const int kstart)
    {
        std::ifstream in(fname, std::ios::in | std::ios::binary);
        int n;
        in.read(reinterpret_cast<char*>(&n),sizeof(n));
        if (n == tmaxi) {
            in.read(reinterpret_cast<char*>(&n),sizeof(n));
            if (n == tmaxj) {
                in.read(reinterpret_cast<char*>(&n),sizeof(n));
                in.read(reinterpret_cast<char*>(&(thedata[kstart*tmaxi*tmaxj])),
                        n*tmaxi*tmaxj*sizeof(double));
                in.close();
            } else {
                std::cerr << "Error reading binary file: tmaxj doesn't match\n";
            }
        } else {
            std::cerr << "Error reading binary file: tmaxi doesn't match\n";
            exit(1);
        }

    }

}
