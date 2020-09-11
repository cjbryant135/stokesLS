#include "uniformmesh2d.h"
#include "numtrait.h"
#include "debug.h"
#include "um2boundary.h"
#include "utility.h"
#include <stdlib.h>

namespace levelset {
        
    UniformMesh2D::UniformMesh2D(STEPSIZE, const int m, const int n,
                                 const int sd, const int si,
                                 const double deltax, const double deltay,
                                 const double z1, const double z2,
                                 UM2_Boundary& thebc)
        : maxi(m), maxj(n), maxk(0), nextw(0), dx(deltax), dy(deltay),
          thedata(NULL), bc(&thebc)
    {
        SetStepSize(m, n, sd, si, deltax, deltay, z1, z2, thebc);
    }

    UniformMesh2D::UniformMesh2D(BOUNDS, const int m, const int n,
                                 const int sd, const int si,
                                 const double xmin, const double xmax,
                                 const double ymin, const double ymax,
                                 UM2_Boundary& thebc)
        : maxi(m), maxj(n), maxk(0), nextw(0), thedata(NULL), bc(&thebc)
    {
        SetBoundsSize(m, n, sd, si, xmin, xmax, ymin, ymax, thebc);
    }

    UniformMesh2D::~UniformMesh2D(void)
    {
        if (thedata) delete[] thedata;
        if (theidata) delete[] theidata;
        if (bc) delete bc;
    }

    void UniformMesh2D::SetStepSize(const int m, const int n,
                                    const int sd, const int si,
                                    const double deltax, const double deltay,
                                    const double zx, const double zy,
                                    UM2_Boundary& thebc)
    {
        maxi = m;
        maxj = n;
        zero[0] = zx;
        zero[1] = zy;
        bc = &thebc;
        bc->SetMesh(this);
        tmaxi = bc->Width(BDRY2D_XLO)+bc->Width(BDRY2D_XHI)+maxi;
        tmaxj = bc->Width(BDRY2D_YLO)+bc->Width(BDRY2D_YHI)+maxj;
        BuildWorkGrid(sd);
        theidata = new int[tmaxi*tmaxj*si];
        for (int i=0; i<tmaxi*tmaxj*si; ++i) theidata[i] = 0;
        dx = deltax;
        dy = deltay;
    }
   

    void UniformMesh2D::SetBoundsSize(const int m, const int n,
                                      const int sd, const int si,
                                      const double xmin, const double xmax,
                                      const double ymin, const double ymax,
                                      UM2_Boundary& thebc)
    {
        maxi = m; 
        maxj = n;
        zero[0] = xmin;
        zero[1] = ymin;
        bc = &thebc;
        bc->SetMesh(this);
        tmaxi = bc->Width(BDRY2D_XLO)+bc->Width(BDRY2D_XHI)+maxi;
        tmaxj = bc->Width(BDRY2D_YLO)+bc->Width(BDRY2D_YHI)+maxj;
        dx = (xmax-xmin)/(m-1);
        dy = (ymax-ymin)/(n-1);
        BuildWorkGrid(sd);
        theidata = new int[tmaxi*tmaxj*si];
        for (int i=0; i<tmaxi*tmaxj*si; ++i) theidata[i] = 0;
    }


    void UniformMesh2D::GetInterface(Interface& seg, const int func,
                                     const double value) 
    {
        enum {FIRST=0, SECOND};
        int chunk = int(sqrt(double(maxi*maxj)));

        if (func > -1) {
            int count = 0;
            seg.Resize(0);
            int ind;
            int i, j;
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j) {
                    idata_(i,j,FIRST) = -1;
                    idata_(i,j,SECOND) = -1;
                }

            for (i=0; i<maxi-1; ++i) 
                for (j=0; j<maxj-1; ++j) {
                    ind =  data_(i,j,func)     >= value ? 1 : 0;
                    ind += data_(i+1,j,func)   >= value ? 2 : 0;
                    ind += data_(i+1,j+1,func) >= value ? 4 : 0;
                    ind += data_(i,j+1,func)   >= value ? 8 : 0;

                    // Make sure seg is large enough to hold new points.

                    switch(ind) {
                    case 0:
                    case 15:
                        break;
                    case 5:
                    case 10:
                        if (count+1 >= seg.Length()) seg.Resize(seg.Length()+chunk);
                        break;
                    default:
                        if (count >= seg.Length()) seg.Resize(seg.Length()+chunk);
                        break;
                    }

                    // Add new points to seg list

                    switch(ind) {
                    case 0: 
                    case 15: // no crossings
                        break;

                    case 1: // SW > value
                        seg[count].SetXY(X(i+frac(value,data_(i,j,func),
                                                  data_(i+1,j,func))), 
                                         Y(j), X(i),
                                         Y(j+frac(value,data_(i,j,func),
                                                  data_(i,j+1,func))),
                                         SFromDown | SToLeft,
                                         seg.datalen);

                        idata_(i,j,FIRST) = count;
                        if (j==0)
                            seg[count].nabor[SLeft] = -1;
                        else
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y2() == Y(j)) {
                                seg[count].nabor[SLeft] = idata_(i,j-1,FIRST);
                                seg[idata_(i,j-1,FIRST)].nabor[SRight] = count;
                            } else {
                                seg[count].nabor[SLeft] = idata_(i,j-1,SECOND);
                                seg[idata_(i,j-1,SECOND)].nabor[SRight] = count;
                            }
                        if (i==0)
                            seg[count].nabor[SRight] = -1;
                        else
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X1() == X(i)) {
                                seg[count].nabor[SRight] = idata_(i-1,j,FIRST);
                                seg[idata_(i-1,j,FIRST)].nabor[SLeft] = count;
                            } else {
                                seg[count].nabor[SRight] = idata_(i-1,j,SECOND);
                                seg[idata_(i-1,j,SECOND)].nabor[SLeft] = count;
                            }

                        ++count;
                        break;

                    case 2: // SE > value
                        seg[count].SetXY(X(i+1),
                                         Y(j+frac(value,data_(i+1,j,func),
                                                  data_(i+1,j+1,func))),
                                         X(i+frac(value,data_(i,j,func),
                                                  data_(i+1,j,func))),
                                         Y(j), SFromRight | SToDown,
                                         seg.datalen);
                                
                        idata_(i,j,FIRST) = count;
                        if (j==0)
                            seg[count].nabor[SRight] = -1;
                        else
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y1() == Y(j)) {
                                seg[count].nabor[SRight] = idata_(i,j-1,FIRST);
                                seg[idata_(i,j-1,FIRST)].nabor[SLeft] = count;
                            } else {
                                seg[count].nabor[SRight] = idata_(i,j-1,SECOND);
                                seg[idata_(i,j-1,SECOND)].nabor[SLeft] = count;
                            }
                        if (i==maxi-2)
                            seg[count].nabor[SLeft] = -1;
                        ++count;
                        break;

                    case 3: // SW,SE > value
                        seg[count].SetXY(X(i+1),
                                         Y(j+frac(value,data_(i+1,j,func),
                                                  data_(i+1,j+1,func))),
                                         X(i),
                                         Y(j+frac(value,data_(i,j,func),
                                                  data_(i,j+1,func))),
                                         SFromRight | SToLeft,
                                         seg.datalen);

                        idata_(i,j,FIRST) = count;
                        if (i==0)
                            seg[count].nabor[SRight] = -1;
                        else
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X1() == X(i)) {
                                seg[count].nabor[SRight] = idata_(i-1,j,FIRST);
                                seg[idata_(i-1,j,FIRST)].nabor[SLeft] = count;
                            } else {
                                seg[count].nabor[SRight] = idata_(i-1,j,SECOND);
                                seg[idata_(i-1,j,SECOND)].nabor[SLeft] = count;
                            }
                        if (i==maxi-2)
                            seg[count].nabor[SLeft] = -1;
                        ++count;
                        break;

                    case 4: // NE > value
                        seg[count].SetXY(X(i+frac(value,data_(i,j+1,func),
                                                  data_(i+1,j+1,func))), 
                                         Y(j+1), X(i+1),
                                         Y(j+frac(value,data_(i+1,j,func),
                                                  data_(i+1,j+1,func))),
                                         SToRight | SFromUp,
                                         seg.datalen);
                                 
                        idata_(i,j,FIRST) = count;
                        if (i==maxi-2)
                            seg[count].nabor[SRight] = -1;
                        if (j==maxj-2)
                            seg[count].nabor[SLeft] = -1;
                        ++count;
                        break;
            
                    case 5: // SW,NE > value
                    {
                        double x[6], y[6];
                        x[0] = X(i+frac(value,data_(i,j,func),data_(i+1,j,func)));
                        y[0] = Y(j);
                        x[1] = X(i+1);
                        y[1] = j+frac(value,data_(i+1,j,func),data_(i+1,j+1,func));
                        x[2] = i+frac(value,data_(i,j+1,func),data_(i+1,j+1,func));
                        y[2] = j+1;
                        x[3] = i;
                        y[3] = j+frac(value,data_(i,j,func),data_(i,j+1,func));
                        int segl,segd;
                        if (j==0) {
                            x[4] = x[0]; y[4] = y[0]; segd = -1;
                        } else {
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y2() == Y(j)) {
                                segd = idata_(i,j-1,FIRST);
                            } else {
                                segd = idata_(i,j-1,SECOND);
                            }
                            x[4] = seg[segd].X2();
                            y[4] = seg[segd].Y2();
                        }
                        if (i==0) {
                            x[5] = x[3]; y[5] = y[3]; segl = -1;
                        } else {
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X1() == X(i)) {
                                segl = idata_(i-1,j,FIRST);
                            } else {
                                segl = idata_(i-1,j,SECOND);
                            }
                            x[5] = seg[segl].X1();
                            y[5] = seg[segl].Y1();
                        }
                        if (fabs((x[0]-x[4])*(y[1]-y[0])-(x[1]-x[0])*(y[0]-y[4]))
                            +fabs((x[3]-x[5])*(y[2]-y[3])-(x[2]-x[3])*(y[3]-y[5]))
                            <fabs((x[0]-x[4])*(y[3]-y[0])-(x[3]-x[0])*(y[0]-y[4]))
                            +fabs((x[3]-x[5])*(y[0]-y[3])-(x[0]-x[3])*(y[3]-y[5]))) {
                            seg[count].SetXY(x[0]*dx+zero[0],y[0]*dy+zero[1],
                                             x[1]*dx+zero[0],y[1]*dy+zero[1],
                                             SToRight | SFromDown,
                                             seg.datalen);
                            idata_(i,j,FIRST) = count;
                            seg[count].nabor[SLeft] = segd;
                            seg[segd].nabor[SRight] = count;
                            if (i==maxi-2)
                                seg[count].nabor[SRight] = -1;
                            ++count;
                            seg[count].SetXY(x[2]*dx+zero[0],y[2]*dy+zero[1],
                                             x[3]*dx+zero[0],y[3]*dy+zero[1],
                                             SToLeft | SFromUp,
                                             seg.datalen);
                            idata_(i,j,SECOND) = count;
                            seg[count].nabor[SRight] = segl;
                            seg[segl].nabor[SLeft] = count;
                            if (j==maxj-2)
                                seg[count].nabor[SLeft] = -1;
                            ++count;
                        } else {
                            seg[count].SetXY(x[0]*dx+zero[0],y[0]*dy+zero[1],
                                             x[3]*dx+zero[0],y[3]*dy+zero[1],
                                             SFromDown | SToLeft,
                                             seg.datalen);
                            idata_(i,j,FIRST) = count;
                            seg[count].nabor[SLeft] = segd;
                            seg[segd].nabor[SRight] = count;
                            seg[count].nabor[SRight] = segl;
                            seg[segl].nabor[SLeft] = count;
                            ++count;
                            seg[count].SetXY(x[2]*dx+zero[0],y[2]*dy+zero[1],
                                             x[1]*dx+zero[0],y[1]*dy+zero[1],
                                             SToRight | SFromUp,
                                             seg.datalen);
                            idata_(i,j,SECOND) = count;
                            if (i==maxi-2)
                                seg[count].nabor[SRight] = -1;
                            if (j==maxj-2)
                                seg[count].nabor[SLeft] = -1;
                            ++count;
                        }
                    }
                    break;
            
                    case 6: // SE,NE > value
                        seg[count].SetXY((i+frac(value,data_(i,j+1,func),
                                                 data_(i+1,j+1,func)))*dx+zero[0],
                                         Y(j+1),
                                         (i+frac(value,data_(i,j,func),
                                                 data_(i+1,j,func)))*dx+zero[0],
                                         Y(j), SFromUp | SToDown,
                                         seg.datalen);
                                
                        idata_(i,j,FIRST) = count;
                        if (j==0)
                            seg[count].nabor[SRight] = -1;
                        else
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y1() == Y(j)) {
                                seg[count].nabor[SRight] = idata_(i,j-1,FIRST);
                                seg[idata_(i,j-1,FIRST)].nabor[SLeft] = count;
                            } else {
                                seg[count].nabor[SRight] = idata_(i,j-1,SECOND);
                                seg[idata_(i,j-1,SECOND)].nabor[SLeft] = count;
                            }
                        if (j==maxj-2)
                            seg[count].nabor[SLeft] = -1;
                        ++count;
                        break;

                    case 7: // SW,SE,NE > value
                        seg[count].SetXY((i+frac(value,data_(i,j+1,func),
                                                 data_(i+1,j+1,func)))*dx+zero[0], 
                                         Y(j+1),
                                         X(i),
                                         (j+frac(value,data_(i,j,func),
                                                 data_(i,j+1,func)))*dy+zero[1],
                                         SToLeft | SFromUp,
                                         seg.datalen);

                        idata_(i,j,FIRST) = count;
                        if (i==0)
                            seg[count].nabor[SRight] = -1;
                        else
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X1() == X(i)) {
                                seg[count].nabor[SRight] = idata_(i-1,j,FIRST);
                                seg[idata_(i-1,j,FIRST)].nabor[SLeft] = count;
                            } else {
                                seg[count].nabor[SRight] = idata_(i-1,j,SECOND);
                                seg[idata_(i-1,j,SECOND)].nabor[SLeft] = count;
                            }
                        if (j==maxj-2)
                            seg[count].nabor[SLeft] = -1;
                        ++count;
                        break;

                    case 8: // NW > value
                        seg[count].SetXY(X(i),
                                         (j+frac(value,data_(i,j,func),
                                                 data_(i,j+1,func)))*dy+zero[1],
                                         (i+frac(value,data_(i,j+1,func),
                                                 data_(i+1,j+1,func)))*dx+zero[0], 
                                         Y(j+1),
                                         SFromLeft | SToUp,
                                         seg.datalen);

                        idata_(i,j,FIRST) = count;
                        if (i==0)
                            seg[count].nabor[SLeft] = -1;
                        else
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X2() == X(i)) {
                                seg[count].nabor[SLeft] = idata_(i-1,j,FIRST);
                                seg[idata_(i-1,j,FIRST)].nabor[SRight] = count;
                            } else {
                                seg[count].nabor[SLeft] = idata_(i-1,j,SECOND);
                                seg[idata_(i-1,j,SECOND)].nabor[SRight] = count;
                            }
                        if (j==maxj-2)
                            seg[count].nabor[SRight] = -1;
                        ++count;
                        break;

                    case 9: // SW,NW > value
                        seg[count].SetXY((i+frac(value,data_(i,j,func),
                                                 data_(i+1,j,func)))*dx+zero[0],
                                         Y(j),
                                         (i+frac(value,data_(i,j+1,func),
                                                 data_(i+1,j+1,func)))*dx+zero[0],
                                         Y(j+1),
                                         SToUp | SFromDown,
                                         seg.datalen);
                                
                        idata_(i,j,FIRST) = count;
                        if (j==0)
                            seg[count].nabor[SLeft] = -1;
                        else
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y2() == Y(j)) {
                                seg[count].nabor[SLeft] = idata_(i,j-1,FIRST);
                                seg[idata_(i,j-1,FIRST)].nabor[SRight] = count;
                            } else {
                                seg[count].nabor[SLeft] = idata_(i,j-1,SECOND);
                                seg[idata_(i,j-1,SECOND)].nabor[SRight] = count;
                            }
                        if (j==maxj-2)
                            seg[count].nabor[SRight] = -1;
                        ++count;
                        break;

                    case 10: // SE,NW > value
                    {   
                        double x[6], y[6];
                        x[0] = (i+frac(value,data_(i,j,func),data_(i+1,j,func)));
                        y[0] = j;
                        x[1] = i+1;
                        y[1] = j+frac(value,data_(i+1,j,func),data_(i+1,j+1,func));
                        x[2] = i+frac(value,data_(i,j+1,func),data_(i+1,j+1,func));
                        y[2] = j+1;
                        x[3] = i;
                        y[3] = j+frac(value,data_(i,j,func),data_(i,j+1,func));
                        int segl,segd;
                        if (j==0) {
                            x[4] = x[0]; y[4] = y[0]; segd = -1;
                        } else {
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y1() == Y(j)) {
                                segd = idata_(i,j-1,FIRST);
                            } else {
                                segd = idata_(i,j-1,SECOND);
                            }
                            x[4] = seg[segd].X1();
                            y[4] = seg[segd].Y1();
                        }
                        if (i==0) {
                            x[5] = x[3]; y[5] = y[3]; segl = -1;
                        } else {
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X2() == X(i)) {
                                segl = idata_(i-1,j,FIRST);
                            } else {
                                segl = idata_(i-1,j,SECOND);
                            }
                            x[5] = seg[segl].X2();
                            y[5] = seg[segl].Y2();
                        }
                        if (fabs((x[0]-x[4])*(y[1]-y[0])-(x[1]-x[0])*(y[0]-y[4]))
                            +fabs((x[3]-x[5])*(y[2]-y[3])-(x[2]-x[3])*(y[3]-y[5]))
                            <fabs((x[0]-x[4])*(y[3]-y[0])-(x[3]-x[0])*(y[0]-y[4]))
                            +fabs((x[3]-x[5])*(y[0]-y[3])-(x[0]-x[3])*(y[3]-y[5]))) {
                            seg[count].SetXY(x[1]*dx+zero[0],y[1]*dy+zero[1],
                                             x[0]*dx+zero[0],y[0]*dy+zero[1],
                                             SFromRight | SToDown,
                                             seg.datalen);
                            idata_(i,j,FIRST) = count;
                            seg[count].nabor[SRight] = segd;
                            seg[segd].nabor[SLeft] = count;
                            if (i==maxi-2)
                                seg[count].nabor[SLeft] = -1;
                            ++count;
                            seg[count].SetXY(x[3]*dx+zero[0],y[3]*dy+zero[1],
                                             x[2]*dx+zero[0],y[2]*dy+zero[1],
                                             SFromLeft | SToUp,
                                             seg.datalen);
                            idata_(i,j,SECOND) = count;
                            seg[count].nabor[SLeft] = segl;
                            seg[segl].nabor[SRight] = count;
                            if (j==maxj-2)
                                seg[count].nabor[SRight] = -1;
                            ++count;
                        } else {
                            seg[count].SetXY(x[3]*dx+zero[0],y[3]*dy+zero[1],
                                             x[0]*dx+zero[0],y[0]*dy+zero[1],
                                             SFromLeft | SToDown,
                                             seg.datalen);
                            idata_(i,j,FIRST) = count;
                            seg[count].nabor[SRight] = segd;
                            seg[segd].nabor[SLeft] = count;
                            seg[count].nabor[SLeft] = segl;
                            seg[segl].nabor[SRight] = count;
                            ++count;
                            seg[count].SetXY(x[1]*dx+zero[0],y[1]*dy+zero[1],
                                             x[2]*dx+zero[0],y[2]*dy+zero[1],
                                             SFromRight | SToUp,
                                             seg.datalen);
                            idata_(i,j,SECOND) = count;
                            if (i==maxi-2)
                                seg[count].nabor[SLeft] = -1;
                            if (j==maxj-2)
                                seg[count].nabor[SRight] = -1;
                            ++count;
                        }
                    }
                    break;
            
                    case 11: // SW,SE,NW > value
                        seg[count].SetXY(X(i+1),
                                         (j+frac(value,data_(i+1,j,func),
                                                 data_(i+1,j+1,func)))*dy+zero[1],
                                         (i+frac(value,data_(i,j+1,func),
                                                 data_(i+1,j+1,func)))*dx+zero[0], 
                                         Y(j+1),
                                         SToUp | SFromRight,
                                         seg.datalen);
                                 
                        idata_(i,j,FIRST) = count;
                        if (i==maxi-2)
                            seg[count].nabor[SLeft] = -1;
                        if (j==maxj-2)
                            seg[count].nabor[SRight] = -1;
                        ++count;
                        break;

                    case 12: // NE,NW > value
                        seg[count].SetXY(X(i),
                                         (j+frac(value,data_(i,j,func),
                                                 data_(i,j+1,func)))*dy+zero[1],
                                         X(i+1),
                                         (j+frac(value,data_(i+1,j,func),
                                                 data_(i+1,j+1,func)))*dy+zero[1],
                                         SFromLeft | SToRight,
                                         seg.datalen);

                        idata_(i,j,FIRST) = count;
                        if (i==0)
                            seg[count].nabor[SLeft] = -1;
                        else
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X2() == X(i)) {
                                seg[count].nabor[SLeft] = idata_(i-1,j,FIRST);
                                seg[idata_(i-1,j,FIRST)].nabor[SRight] = count;
                            } else {
                                seg[count].nabor[SLeft] = idata_(i-1,j,SECOND);
                                seg[idata_(i-1,j,SECOND)].nabor[SRight] = count;
                            }
                        if (i==maxi-2)
                            seg[count].nabor[SRight] = -1;
                        ++count;
                        break;

                    case 13: // SW,NE,NW > value
                        seg[count].SetXY((i+frac(value,data_(i,j,func),
                                                 data_(i+1,j,func)))*dx+zero[0],
                                         Y(j),
                                         X(i+1),
                                         (j+frac(value,data_(i+1,j,func),
                                                 data_(i+1,j+1,func)))*dy+zero[1],
                                         SToRight | SFromDown,
                                         seg.datalen);
                                
                        idata_(i,j,FIRST) = count;
                        if (j==0)
                            seg[count].nabor[SLeft] = -1;
                        else
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y2() == Y(j)) {
                                seg[count].nabor[SLeft] = idata_(i,j-1,FIRST);
                                seg[idata_(i,j-1,FIRST)].nabor[SRight] = count;
                            } else {
                                seg[count].nabor[SLeft] = idata_(i,j-1,SECOND);
                                seg[idata_(i,j-1,SECOND)].nabor[SRight] = count;
                            }
                        if (i==maxi-2)
                            seg[count].nabor[SRight] = -1;
                        ++count;
                        break;

                    case 14: // SE,NE,NW > value
                        seg[count].SetXY(X(i),
                                         (j+frac(value,data_(i,j,func),
                                                 data_(i,j+1,func)))*dy+zero[1],
                                         (i+frac(value,data_(i,j,func),
                                                 data_(i+1,j,func)))*dx+zero[0], 
                                         Y(j),
                                         SFromLeft | SToDown,
                                         seg.datalen);

                        idata_(i,j,FIRST) = count;
                        if (j==0)
                            seg[count].nabor[SRight] = -1;
                        else
                            if (idata_(i,j-1,SECOND) == -1 || 
                                seg[idata_(i,j-1,FIRST)].Y1() == Y(j)) {
                                seg[count].nabor[SRight] = idata_(i,j-1,FIRST);
                                seg[idata_(i,j-1,FIRST)].nabor[SLeft] = count;
                            } else {
                                seg[count].nabor[SRight] = idata_(i,j-1,SECOND);
                                seg[idata_(i,j-1,SECOND)].nabor[SLeft] = count;
                            }
                        if (i==0)
                            seg[count].nabor[SLeft] = -1;
                        else
                            if (idata_(i-1,j,SECOND) == -1 || 
                                seg[idata_(i-1,j,FIRST)].X2() == X(i)) {
                                seg[count].nabor[SLeft] = idata_(i-1,j,FIRST);
                                seg[idata_(i-1,j,FIRST)].nabor[SRight] = count;
                            } else {
                                seg[count].nabor[SLeft] = idata_(i-1,j,SECOND);
                                seg[idata_(i-1,j,SECOND)].nabor[SRight] = count;
                            }

                        ++count;
                        break;
                    }
                }

            seg.Resize(count);
        } else {
            seg.Resize(0);
        }
    }


    void UniformMesh2D::Interp(Interface& seg, const int res, const int k)
    {
        int len = seg.Length();
        double xloc, yloc;
        int i;
        for (i=0; i<len; ++i) {
            seg[i].MidPoint(xloc, yloc);
            seg[i][res] = Interp(xloc, yloc, k);
        }
    }

    double UniformMesh2D::Interp2(const double& x, const double& y, const int k) const
    {
        int i, j;
        double ifrac, jfrac;
        LocToIndex(x, y, i, j, ifrac, jfrac);
        return Interp2(i, j, k, ifrac, jfrac);
    }

    double UniformMesh2D::Interp2(const int i, const int j, const int k, const double& ifrac,
                                  const double& jfrac) const
    {
        Bicubic p;
        p.BuildwDeriv(data_(i-1,j-1,k),data_(i,j-1,k),data_(i+1,j-1,k),
                      data_(i+2,j-1,k),data_(i-1,j,k),data_(i,j,k),data_(i+1,j,k),
                      data_(i+2,j,k),data_(i-1,j+1,k),data_(i,j+1,k),
                      data_(i+1,j+1,k),data_(i+2,j+1,k),data_(i-1,j+2,k),
                      data_(i,j+2,k),data_(i+1,j+2,k),data_(i+2,j+2,k),dx,dy);
        return p(ifrac*dx, jfrac*dy);
    }

    Bicubic UniformMesh2D::Interp2(const int i, const int j, const int k) const
    {
        Bicubic p;
        p.BuildwDeriv(data_(i-1,j-1,k),data_(i,j-1,k),data_(i+1,j-1,k),
                      data_(i+2,j-1,k),data_(i-1,j,k),data_(i,j,k),data_(i+1,j,k),
                      data_(i+2,j,k),data_(i-1,j+1,k),data_(i,j+1,k),
                      data_(i+1,j+1,k),data_(i+2,j+1,k),data_(i-1,j+2,k),
                      data_(i,j+2,k),data_(i+1,j+2,k),data_(i+2,j+2,k),dx,dy);
        return p;
    }


    void UniformMesh2D::LocToIndex(const double x, const double y, int& i, int& j, 
                                   double& ifrac, double& jfrac, const char round) const
    {
        if (round) {
            i = floor((x-zero[0])/dx+0.5);
            j = floor((y-zero[1])/dy+0.5);
            ifrac = (x-zero[0]-i*dx)/dx;
            jfrac = (y-zero[1]-j*dy)/dy;
        } else {
            i = floor((x-zero[0])/dx);
            j = floor((y-zero[1])/dy);
            ifrac = (x-zero[0]-i*dx)/dx;
            jfrac = (y-zero[1]-j*dy)/dy;
        }
    }

    void UniformMesh2D::LocToIndex(const double x, const double y, int& i, int& j, 
                                   const char round) const
    {
        if (round) {
            i = int((x-zero[0])/dx+0.5);
            j = int((y-zero[1])/dy+0.5);
        } else {
            i = int((x-zero[0])/dx);
            j = int((y-zero[1])/dy);
        }
    }


    void UniformMesh2D::BuildWorkGrid(int num)
    {
        thedata = new double[tmaxi*tmaxj*num];
        datastart = Index(bc->Width(BDRY2D_XLO), bc->Width(BDRY2D_YLO));
        maxk = num;
        for (int i=0; i<tmaxi*tmaxj*num; ++i) thedata[i] = 0.;
    }

    void UniformMesh2D::CopyWorkGrid(int to, int from)
    {
        double* datato = &(thedata[Index(0,0,to)]);
        double* datafrom = &(thedata[Index(0,0,from)]);
        memcpy(datato, datafrom, tmaxi*tmaxj*sizeof(double));
    }

    void UniformMesh2D::CopyIntGrid(int to, int from)
    {
        int* datato = &(theidata[Index(0,0,to)]);
        int* datafrom = &(theidata[Index(0,0,from)]);
        memcpy(datato, datafrom, tmaxi*tmaxj*sizeof(int));
    }

    inline double triangle_area(double x1, double y1, double x2, double y2,
                                double x3, double y3)
    {
        return fabs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2.;
    }

    double UniformMesh2D::Area(const int k, const double val) 
    {
        double area = 0.;
   
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
                int d = data_(i,j,k) >= val ? 1 : 0;
                d += data_(i+1,j,k) >= val ? 2 : 0;
                d += data_(i+1,j+1,k) >= val ? 4 : 0;
                d += data_(i,j+1,k) >= val ? 8 : 0;

                switch (d) {
                case 0:
                    break;
                case 15:
                    area += dx*dy;
                    break;
                case 1:
                    area += dx*frac(val,data_(i,j,k),data_(i+1,j,k))
                        *dy*frac(val,data_(i,j,k),data_(i,j+1,k))/2.;
                    break;
                case 2:
                    area +=  dx*frac(val,data_(i+1,j,k),data_(i,j,k))
                        *dy*frac(val,data_(i+1,j,k),data_(i+1,j+1,k))/2.;
                    break;
                case 3:
                    area += dx*dy*(frac(val,data_(i,j,k),data_(i,j+1,k))
                                   +frac(val,data_(i+1,j,k),data_(i+1,j+1,k)))/2.;
                    break;
                case 4:
                    area += dx*frac(val,data_(i+1,j+1,k),data_(i,j+1,k))
                        *dy*frac(val,data_(i+1,j+1,k),data_(i+1,j,k))/2.;
                    break;
                case 5:
                {
                    double x1 = frac(val,data_(i,j,k),data_(i+1,j,k));
                    double y1 = frac(val,data_(i,j,k),data_(i,j+1,k));
                    double x2 = frac(val,data_(i,j+1,k),data_(i+1,j+1,k));
                    double y2 = frac(val,data_(i+1,j,k),data_(i+1,j+1,k));
                    double cx = (-y1*(x2-x1)-x1)/((y2-y1)*(x2-x1)-1);
                    double cy = (-x1*(y2-y1)-y1)/((y2-y1)*(x2-x1)-1);
                    area += triangle_area(0.,0.,x1*dx,0.,cx*dx,cy*dy)
                        +triangle_area(0.,0.,cx*dx,cy*dy,0,y1*dy)
                        +triangle_area(cx*dx,cy*dy,dx,y2*dy,dx,dy)
                        +triangle_area(cx*dx,cy*dy,dx,dy,x2*dx,dy);
                }
                break;
                case 6:
                    area += dx*dy*(frac(val,data_(i+1,j,k),data_(i,j,k))
                                   +frac(val,data_(i+1,j+1,k),data_(i,j+1,k)))/2.;
                    break;
                case 7:
                    area += dx*dy-dx*frac(val,data_(i,j+1,k),data_(i+1,j+1,k))
                        *dy*frac(val,data_(i,j+1,k),data_(i,j,k))/2.;
                    break;
                case 8:
                    area += dx*frac(val,data_(i,j+1,k),data_(i+1,j+1,k))
                        *dy*frac(val,data_(i,j+1,k),data_(i,j,k))/2.;
                    break;
                case 9:
                    area += dx*dy*(frac(val,data_(i,j,k),data_(i+1,j,k))
                                   +frac(val,data_(i,j+1,k),data_(i+1,j+1,k)))/2.;
                    break;
                case 10:
                {
                    double x1 = frac(val,data_(i,j,k),data_(i+1,j,k));
                    double y1 = frac(val,data_(i,j,k),data_(i,j+1,k));
                    double x2 = frac(val,data_(i,j+1,k),data_(i+1,j+1,k));
                    double y2 = frac(val,data_(i+1,j,k),data_(i+1,j+1,k));
                    double cx = (-y1*(x2-x1)-x1)/((y2-y1)*(x2-x1)-1);
                    double cy = (-x1*(y2-y1)-y1)/((y2-y1)*(x2-x1)-1);
                    area += dx*dy-triangle_area(0.,0.,x1*dx,0.,cx*dx,cy*dy)
                        -triangle_area(0.,0.,cx*dx,cy*dy,0,y1*dy)
                        -triangle_area(cx*dx,cy*dy,dx,y2*dy,dx,dy)
                        -triangle_area(cx*dx,cy*dy,dx,dy,x2*dx,dy);
                }
                break;
                case 11:
                    area += dx*dy-dx*frac(val,data_(i+1,j+1,k),data_(i,j+1,k))
                        *dy*frac(val,data_(i+1,j+1,k),data_(i+1,j,k))/2.;
                    break;
                case 12:
                    area += dx*dy*(frac(val,data_(i,j+1,k),data_(i,j,k))
                                   +frac(val,data_(i+1,j+1,k),data_(i+1,j,k)))/2.;
                    break;
                case 13:
                    area += dx*dy-dx*frac(val,data_(i+1,j,k),data_(i,j,k))
                        *dy*frac(val,data_(i+1,j,k),data_(i+1,j+1,k))/2.;
                    break;
                case 14:
                    area += dx*dy-dx*frac(val,data_(i,j,k),data_(i+1,j,k))
                        *dy*frac(val,data_(i,j,k),data_(i,j+1,k))/2.;
                    break;
                }
            }
        return area;
    }

    double UniformMesh2D::IntegrateTrap(const int k) 
    {
        double sum = 0.;
   
        sum = data_(0,0,k)/4.;
        for (int j=1; j<maxj-1; ++j) 
            sum += data_(0,j,k)/2.;
        sum += data_(0,maxj-1,k)/4.;
        for (int i=1; i<maxi-1; ++i) {
            sum += data_(i,0,k)/2.; 
            for (int j=1; j<maxj-1; ++j)
                sum += data_(i,j,k);
            sum += data_(i,maxj-1,k)/2.;
        }
        sum += data_(maxi-1,0,k)/4.;
        for (int j=1; j<maxj-1; ++j)
            sum += data_(maxi-1,j,k)/2.;
        sum += data_(maxi-1,maxj-1,k)/4.;

        return sum*dx*dy;
    }

    double UniformMesh2D::IntegrateBicubic(const int k) 
    {
        double sum = 0.;
        Bicubic p;
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
                p.BuildwDeriv(data_(i-i,j-1,k),data_(i,j-1,k),data_(i+1,j-1,k),data_(i+2,j-1,k),
                              data_(i-i,j,k),data_(i,j,k),data_(i+1,j,k),data_(i+2,j,k),
                              data_(i-i,j+1,k),data_(i,j+1,k),data_(i+1,j+1,k),data_(i+2,j+1,k),
                              data_(i-i,j+2,k),data_(i,j+2,k),data_(i+1,j+2,k),data_(i+2,j+2,k),
                              dx, dy);
                sum += p.Integral();
            }

        return sum;
    }

#if 0
    void UniformMesh2D::CenterOfMass(const int k, const double val,
                                     double& cx, double& cy) 
    {
        double area = 0.;
        double e[4];
        cx = 0.;
        cy = 0.;
   
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j) {
                int d = data_(i,j,k) >= val ? 1 : 0;
                d += data_(i+1,j,k) >= val ? 2 : 0;
                d += data_(i+1,j+1,k) >= val ? 4 : 0;
                d += data_(i,j+1,k) >= val ? 8 : 0;

                switch (d) {
                case 0:
                    break;
                case 15:
                    area += dx*dy;
                    cx += dx*dy/2.*(dx+2.*X(i));
                    cy += dx*dy/2.*(dy+2.*Y(j));
                    break;
                case 1:
                    e[0] = dx*frac(val,data_(i,j,k),data_(i+1,j,k));
                    e[1] = dy*frac(val,data_(i,j,k),data_(i,j+1,k));
                    area += e[0]*e[1]/2.;
                    cx += e[0]*e[1]*(e[0]+3*X(i))/6.;
                    cy += e[0]*e[1]*(e[1]+3*Y(j))/6.;
                    break;
                case 2:
                    e[0] = dx*frac(val,data_(i+1,j,k),data_(i,j,k));
                    e[1] = dy*frac(val,data_(i+1,j,k),data_(i+1,j+1,k));
                    area +=  e[0]*e[1]/2.;
                    cx += e[0]*e[1]*(3*X(i+1)-e[0])/6.;
                    cy += e[0]*e[1]*(e[1]+3*Y(j))/6.;
                    break;
                case 3:
                    e[0] = dy*frac(val,data_(i,j,k),data_(i,j+1,k));
                    e[1] = dy*frac(val,data_(i+1,j,k),data_(i+1,j+1,k));
                    area += dx*(e[0]+e[1])/2.;
                    cx += dx*(3*X(i)*(e[0]+e[1])+dx*(e[0]+2*e[1]))/6.;
                    cy += dx*(3*Y(j)*(e[0]+e[1])+sqr(e[0]+e[1])-e[0]*e[1])/6.;
                    break;
                case 4:
                    area += dx*frac(val,data_(i+1,j+1,k),data_(i,j+1,k))
                        *dy*frac(val,data_(i+1,j+1,k),data_(i+1,j,k))/2.;
                    break;
                case 5:
                {
                    double x1 = frac(val,data_(i,j,k),data_(i+1,j,k));
                    double y1 = frac(val,data_(i,j,k),data_(i,j+1,k));
                    double x2 = frac(val,data_(i,j+1,k),data_(i+1,j+1,k));
                    double y2 = frac(val,data_(i+1,j,k),data_(i+1,j+1,k));
                    double cx = (-y1*(x2-x1)-x1)/((y2-y1)*(x2-x1)-1);
                    double cy = (-x1*(y2-y1)-y1)/((y2-y1)*(x2-x1)-1);
                    area += triangle_area(0.,0.,x1*dx,0.,cx*dx,cy*dy)
                        +triangle_area(0.,0.,cx*dx,cy*dy,0,y1*dy)
                        +triangle_area(cx*dx,cy*dy,dx,y2*dy,dx,dy)
                        +triangle_area(cx*dx,cy*dy,dx,dy,x2*dx,dy);
                }
                break;
                case 6:
                    area += dx*dy*(frac(val,data_(i+1,j,k),data_(i,j,k))
                                   +frac(val,data_(i+1,j+1,k),data_(i,j+1,k)))/2.;
                    break;
                case 7:
                    area += dx*dy-dx*frac(val,data_(i,j+1,k),data_(i+1,j+1,k))
                        *dy*frac(val,data_(i,j+1,k),data_(i,j,k))/2.;
                    break;
                case 8:
                    e[0] = dx*frac(val,data_(i,j+1,k),data_(i+1,j+1,k));
                    e[1] = dy*frac(val,data_(i,j+1,k),data_(i,j,k));
                    area += e[0]*e[1]/2.;
                    cx += e[0]*e[1]*(3*X(i)+e[0])/6.;
                    cy += e[0]*e[1]*(3*Y(j+1)-e[1])/6.;
                    break;
                case 9:
                    e[0] = dx*frac(val,data_(i,j,k),data_(i+1,j,k));
                    e[1] = dx*frac(val,data_(i,j+1,k),data_(i+1,j+1,k));
                    area += dy*(e[0]+e[1])/2.;
                    cx += dy*(3*X(i)*(e[0]+e[1])+sqr(e[0]+e[1])-e[0]*e[1])/6.;
                    cy += dy*(3*Y(j)*(e[0]+e[1])+dy*(e[0]+2*e[1]))/6.;
                    break;
                case 10:
                    e[0] = dx*frac(val,data_(i+1,j,k),data_(i,j,k));
                    e[1] = dy*frac(val,data_(i+1,j,k),data_(i+1,j+1,k));
                    e[2] = dx*frac(val,data_(i,j+1,k),data_(i+1,j+1,k));
                    e[3] = dy*frac(val,data_(i,j+1,k),data_(i,j,k));
                    if (e[0]+e[1]+e[2]+e[3] > dx+dy) {
                        area += ((e[1]+e[3])*dx+(e[0]+e[2])*dy-e[1]*e[2]-e[0]*e[3])/2.;
                        cx += X(i)*(dx*(e[1]+e[3])+dy*(e[0]+e[2])
                                    -e[0]*e[3]-e[1]*e[2])/2.+dx*dx*(2*e[1]+e[3])/6.
                            +dx*dy*(2*e[0]+e[2])/6.-dx*(2*e[0]*e[3]-e[1]*e[2])/6.
                            +dy*(e[2]*e[2]-e[0]*e[0])/6.
                            +(e[0]*e[0]*e[3]-e[1]*e[2]*e[2])/6.;
                        cy += Y(j)*(dx*(e[1]+e[3])+dy*(e[0]+e[2])
                                    -e[1]*e[2]-e[0]*e[3])/2.+dy*dy*(2*e[2]+e[0])/6.
                            +dx*dy*(2*e[3]+e[1])/6.-dy*(2*e[0]*e[3]+e[1]*e[2])/6.
                            +dx*(e[1]*e[1]-e[3]*e[3])/6.
                            +(e[0]*e[3]*e[3]-e[1]*e[1]*e[2])/6.;
                    }
                    else {
                        area += (e[0]*e[1]+e[2]*e[3])/2.;
                        cx += e[2]*e[3]*(3*X(i)+e[2])/6.+e[0]*e[1]*(3*X(i+1)-e[0])/6.;
                        cy += e[2]*e[3]*(3*Y(j+1)-e[3])/6.+e[0]*e[1]*(e[1]+3*Y(j))/6.;
                    }
                    break;
                case 11:
                    e[0] = dx*frac(val,data_(i+1,j+1,k),data_(i,j+1,k));
                    e[1] = dy*frac(val,data_(i+1,j+1,k),data_(i+1,j,k));
                    area += dx*dy-e[0]*e[1]/2.;
                    cx += X(i)*(dx*dy-e[0]*e[1]/2.)
                        break;
                case 12:
                    area += dx*dy*(frac(val,data_(i,j+1,k),data_(i,j,k))
                                   +frac(val,data_(i+1,j+1,k),data_(i+1,j,k)))/2.;
                    break;
                case 13:
                    area += dx*dy-dx*frac(val,data_(i+1,j,k),data_(i,j,k))
                        *dy*frac(val,data_(i+1,j,k),data_(i+1,j+1,k))/2.;
                    break;
                case 14:
                    area += dx*dy-dx*frac(val,data_(i,j,k),data_(i+1,j,k))
                        *dy*frac(val,data_(i,j,k),data_(i,j+1,k))/2.;
                    break;
                }
            }
        return area;
    }
#endif

    char UniformMesh2D::IsCrossing(const double value, const double y0, const double y1, double& f)
    {
        if (value-y0 >= 0 && value-y0 <= y1-y0) {
            f = frac(value, y0, y1);
            return true;
        } else
            return false;
    }

            
    char UniformMesh2D::TestValidity(const char* file, const int line,
                                     const int kind) const
    {
        char ans = true;
   
        if (kind >= 0) {
            for (int iii=0; iii<maxi; ++iii)
                for (int jjj=0; jjj<maxj; ++jjj) {
                    if (!finite(data_(iii,jjj,kind))) {
                        std::cerr << "Warning: Invalid data at (" << iii << ',' << jjj
                                  << ',' << "kFunc) in " << file
                                  << ':' << line << '\n';
                        ans = false;
                    }
                }
        } else {
            for (int k=0; k<maxk; ++k) 
                for (int iii=0; iii<maxi; ++iii)
                    for (int jjj=0; jjj<maxj; ++jjj) {
                        if (!finite(data_(iii,jjj,k))) {
                            std::cerr << "Warning: Invalid data at (" << iii
                                      << ',' << jjj << ',' << k
                                      << ") in " << file << ':' << line << '\n';
                            ans = false;
                        }
                    }
        }
        if (ans) std::cerr << "Data OK in " << file << ':' << line << '\n';
        return ans;
    }

    void UniformMesh2D::MinMax(const int k, double& themin, double& themax) const
    {
        themin = data_(0,0,k);
        themax = data_(0,0,k);
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j) {
                if (data_(i,j,k) < themin) themin = data_(i,j,k);
                if (data_(i,j,k) > themax) themax = data_(i,j,k); 
            }
    }

}
