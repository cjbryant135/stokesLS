#include "uniformmesh3d.h"
#include "numtrait.h"
#include "debug.h"
#include "um3boundary.h"
#include "interface3d.h"
#include <stdlib.h>

namespace levelset {
        
    UniformMesh3D::UniformMesh3D(STEPSIZE, const int m, const int n, const int o,
                                 const int sd, const int si,
                                 const double deltax, const double deltay,
                                 const double deltaz, const double z1,
                                 const double z2, const double z3,
                                 UM3_Boundary& thebc)
        : maxi(m), maxj(n), maxk(o), maxl(0), dx(deltax), dy(deltay), dz(deltaz),
          thedata(NULL), bc(&thebc)
    {
        SetStepSize(m, n, o, sd, si, deltax, deltay, deltaz, z1, z2, z3, thebc);
    }

    UniformMesh3D::UniformMesh3D(BOUNDS, const int m, const int n, const int o,
                                 const int sd, const int si,
                                 const double xmin, const double xmax,
                                 const double ymin, const double ymax,
                                 const double zmin, const double zmax,
                                 UM3_Boundary& thebc)
        : maxi(m), maxj(n), maxk(o), maxl(0), thedata(NULL), bc(&thebc)
    {
        SetBoundsSize(m, n, o, sd, si, xmin, xmax, ymin, ymax, zmin, zmax, thebc);
    }

    UniformMesh3D::~UniformMesh3D(void)
    {
        if (thedata) delete[] thedata;
        if (theidata) delete[] theidata;
    }

    void UniformMesh3D::SetStepSize(const int m, const int n, const int o,
                                    const int sd, const int si,
                                    const double deltax, const double deltay,
                                    const double deltaz, const double zx,
                                    const double zy, const double zz,
                                    UM3_Boundary& thebc)
    {
        maxi = m;
        maxj = n;
        maxk = o;
        zero[0] = zx;
        zero[1] = zy;
        zero[2] = zz;
        bc = &thebc;
        bc->SetMesh(this);
        tmaxi = bc->Width(BDRY3D_XLO)+bc->Width(BDRY3D_XHI)+maxi;
        tmaxj = bc->Width(BDRY3D_YLO)+bc->Width(BDRY3D_YHI)+maxj;
        tmaxk = bc->Width(BDRY3D_ZLO)+bc->Width(BDRY3D_ZHI)+maxk;
        BuildWorkGrid(sd);
        theidata = new int[tmaxi*tmaxj*tmaxk*si];
        dx = deltax;
        dy = deltay;
        dz = deltaz;
    }
   

    void UniformMesh3D::SetBoundsSize(const int m, const int n, const int o,
                                      const int sd, const int si,
                                      const double xmin, const double xmax,
                                      const double ymin, const double ymax,
                                      const double zmin, const double zmax,
                                      UM3_Boundary& thebc)
    {
        maxi = m; 
        maxj = n;
        maxk = o;
        zero[0] = xmin;
        zero[1] = ymin;
        zero[2] = zmin;
        bc = &thebc;
        bc->SetMesh(this);
        tmaxi = bc->Width(BDRY3D_XLO)+bc->Width(BDRY3D_XHI)+maxi;
        tmaxj = bc->Width(BDRY3D_YLO)+bc->Width(BDRY3D_YHI)+maxj;
        tmaxk = bc->Width(BDRY3D_ZLO)+bc->Width(BDRY3D_ZHI)+maxk;
        dx = (xmax-xmin)/(m-1);
        dy = (ymax-ymin)/(n-1);
        dz = (zmax-zmin)/(o-1);
        BuildWorkGrid(sd);
        theidata = new int[tmaxi*tmaxj*tmaxk*si];
    }

    enum {
        LYLZ=0x0001, HYLZ=0x0002, HYHZ=0x0004, LYHZ=0x0008, LXLZ=0x0010,
        HXLZ=0x0020, HXHZ=0x0040, LXHZ=0x0080, LXLY=0x0100, HXLY=0x0200,
        HXHY=0x0400, LXHY=0x0800
    };

    inline int over_value(const double grid_value, const double levelset) 
    {return ((grid_value > levelset) ? 1 : 0);}


    void UniformMesh3D::GetInterface(Interface3D& surf, const int func,
                                     const int itmp1, const int itmp2, 
                                     const int itmp3, const double value) 
    {
        int i, j, k;
        double f;
        int chunk = maxi;
        surf.ncount = 0;
        
        // Find all nodes
        
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    idata_(i,j,k,itmp1) = -1;
                    idata_(i,j,k,itmp2) = -1;
                    idata_(i,j,k,itmp3) = -1;
                }
        
        for (i=0; i<maxi-1; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) 
                    if (IsCrossing(value, data_(i,j,k,func), data_(i+1,j,k,func), f)) {
                        if (surf.ncount+1 >= surf.node.Length())
                            surf.node.Resize(surf.node.Length()+chunk);
                        surf.node[surf.ncount].Set(X(i+f),Y(j),Z(k));
                        idata_(i,j,k,itmp1) = surf.ncount++;
                    }

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj-1; ++j)
                for (k=0; k<maxk; ++k) 
                    if (IsCrossing(value, data_(i,j,k,func), data_(i,j+1,k,func), f)) {
                        if (surf.ncount+1 >= surf.node.Length())
                            surf.node.Resize(surf.node.Length()+chunk);
                        surf.node[surf.ncount].Set(X(i),Y(j+f),Z(k));
                        idata_(i,j,k,itmp2) = surf.ncount++;
                    }

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk-1; ++k) 
                    if (IsCrossing(value, data_(i,j,k,func), data_(i,j,k+1,func), f)) {
                        if (surf.ncount+1 >= surf.node.Length())
                            surf.node.Resize(surf.node.Length()+chunk);
                        surf.node[surf.ncount].Set(X(i),Y(j),Z(k+f));
                        idata_(i,j,k,itmp3) = surf.ncount++;
                    }
                                
        // Connect the nodes
        
        int box_value;
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

        surf.ecount = 0;
        chunk = surf.ncount;
        int nodeid[3];
        
        for (i=0; i<maxi-1; ++i) 
            for (j=0; j<maxj-1; ++j) 
                for (k=0; k<maxk-1; ++k) {
               
                    box_value = over_value(data_(i+1,j+1,k+1,func), value);
                    box_value = 2*box_value +over_value(data_(i+1,j+1,k,func), value);
                    box_value = 2*box_value +over_value(data_(i+1,j,k+1,func),value);
                    box_value = 2*box_value +over_value(data_(i+1,j,k,func), value);
                    box_value = 2*box_value +over_value(data_(i,j+1,k+1,func), value);
                    box_value = 2*box_value +over_value(data_(i,j+1,k,func), value);
                    box_value = 2*box_value +over_value(data_(i,j,k+1,func),value);
                    box_value = 2*box_value +over_value(data_(i,j,k,func), value);
                    if (box_value > 127) box_value = 255 - box_value;

                    if (box_value != 0) {
                        for (int m=1; m<=lookup[box_value-1][0]; ++m) {
                                        
                            if (surf.ecount+5 >= surf.elem.Length())
                                surf.elem.Resize(surf.elem.Length()+chunk);
                                        
                            int n;
                            n = 0;
                            if (lookup[box_value-1][m] & LYLZ)
                                nodeid[n++] = idata_(i,j,k,itmp1);
                            if (lookup[box_value-1][m] & HYLZ) 
                                nodeid[n++] = idata_(i,j+1,k,itmp1);
                            if (lookup[box_value-1][m] & HYHZ) 
                                nodeid[n++] = idata_(i,j+1,k+1,itmp1);
                            if (lookup[box_value-1][m] & LYHZ) 
                                nodeid[n++] = idata_(i,j,k+1,itmp1);
                            if (lookup[box_value-1][m] & LXLZ)
                                nodeid[n++] = idata_(i,j,k,itmp2);
                            if (lookup[box_value-1][m] & HXLZ) 
                                nodeid[n++] = idata_(i+1,j,k,itmp2);
                            if (lookup[box_value-1][m] & HXHZ) 
                                nodeid[n++] = idata_(i+1,j,k+1,itmp2);
                            if (lookup[box_value-1][m] & LXHZ) 
                                nodeid[n++] = idata_(i,j,k+1,itmp2);
                            if (lookup[box_value-1][m] & LXLY) 
                                nodeid[n++] = idata_(i,j,k,itmp3);
                            if (lookup[box_value-1][m] & HXLY) 
                                nodeid[n++] = idata_(i+1,j,k,itmp3);
                            if (lookup[box_value-1][m] & HXHY) 
                                nodeid[n++] = idata_(i+1,j+1,k,itmp3);
                            if (lookup[box_value-1][m] & LXHY) 
                                nodeid[n++] = idata_(i,j+1,k,itmp3);
                            if (n == 3) {
                                surf.elem[surf.ecount++].Set(nodeid);
                            } else {
                                std::cout << "ERROR: faulty data table for "
                                          << box_value << '\n';
                            };
                        }
                    }
                }
    }


    void UniformMesh3D::LocToIndex(const double x, const double y, const double z,
                                   int& i, int& j, int& k,
                                   double& ifrac, double& jfrac, double& kfrac,
                                   const char round) const
    {
        if (round) {
            i = int((x-zero[0])/dx+0.5);
            j = int((y-zero[1])/dy+0.5);
            k = int((z-zero[2])/dz+0.5);
            ifrac = (x-zero[0]-i*dx)/dx;
            jfrac = (y-zero[1]-j*dy)/dy;
            kfrac = (z-zero[2]-k*dz)/dz;
        } else {
            i = int((x-zero[0])/dx);
            j = int((y-zero[1])/dy);
            k = int((z-zero[2])/dz);
            ifrac = (x-zero[0]-i*dx)/dx;
            jfrac = (y-zero[1]-j*dy)/dy;
            kfrac = (z-zero[2]-k*dz)/dz;
        }
    }

    void UniformMesh3D::LocToIndex(const double x, const double y, const double z,
                                   int& i, int& j, int& k, const char round) const
    {
        if (round) {
            i = int((x-zero[0])/dx+0.5);
            j = int((y-zero[1])/dy+0.5);
            k = int((z-zero[2])/dz+0.5);
        } else {
            i = int((x-zero[0])/dx);
            j = int((y-zero[1])/dy);
            k = int((z-zero[2])/dz);
        }
    }

    double UniformMesh3D::Interp2(const double x, const double y, const double z, const int l) const
    {
        int i, j, k;
        double ifrac, jfrac, kfrac;
        LocToIndex(x, y, z, i, j, k, ifrac, jfrac, kfrac);
        return Interp2(i, j, k, l, ifrac, jfrac, kfrac);
    }

    double UniformMesh3D::Interp2(const int i, const int j, const int k, const int l, const double& ifrac,
                                  const double& jfrac, const double& kfrac) const
    {
        Tricubic p;
        double f[4][4][4];
   
        for (int ii=-1; ii<3; ++ii)
            for (int jj=-1; jj<3; ++jj)
                for (int kk=-1; kk<3; ++kk)
                    f[ii+1][jj+1][kk+1] = data_(i+ii,j+jj,k+kk,l);
                        
        p.BuildwDeriv(f,dx,dy,dz);
        return p(ifrac*dx, jfrac*dy, kfrac*dz);
    }

    Tricubic UniformMesh3D::Interp2(const int i, const int j, const int k, const int l) const
    {
        Tricubic p;
        double f[4][4][4];
   
        for (int ii=-1; ii<3; ++ii)
            for (int jj=-1; jj<3; ++jj)
                for (int kk=-1; kk<3; ++kk)
                    f[ii+1][jj+1][kk+1] = data_(i+ii,j+jj,k+kk,l);
                        
        p.BuildwDeriv(f,dx,dy,dz);
        return p;
    }



    void UniformMesh3D::BuildWorkGrid(int num)
    {
        thedata = new double[tmaxi*tmaxj*tmaxk*num];
        datastart = Index(bc->Width(BDRY3D_XLO), bc->Width(BDRY3D_YLO),
                          bc->Width(BDRY3D_ZLO));
        maxl = num;
    }

    void UniformMesh3D::CopyWorkGrid(int to, int from)
    {
        double* datato = &(thedata[Index(0,0,0,to)]);
        double* datafrom = &(thedata[Index(0,0,0,from)]);
        memcpy(datato, datafrom, tmaxi*tmaxj*tmaxk*sizeof(double));
    }

    inline double triangle_area(double x1, double y1, double x2, double y2,
                                double x3, double y3)
    {
        return fabs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2.;
    }

    double CalcRot_(const double u000, const double u100, const double u110,
                    const double u010, const double u001, const double u101,
                    const double u111, const double u011, const double lvl,
                    int r)
    {
        double  deltav;
        double  a[4], b[4], c[4];

        switch(r) {
        case 1:         /*      empty/full box  */
            deltav = 0.;
            break;
        
        case 2:
            a[0] = frac(lvl,u000,u100);
            b[0] = frac(lvl,u000,u010);
            c[0] = frac(lvl,u000,u001);
            deltav = a[0]*b[0]*c[0]/6.;
            break;
        
        case 3:
            a[0] = frac(lvl,u000,u100);
            a[1] = frac(lvl,u001,u101);
            b[0] = frac(lvl,u000,u010);
            b[1] = frac(lvl,u001,u011);
            deltav = (a[0]*b[1]+2.*a[0]*b[0]+2.*a[1]*b[1]+a[1]*b[0])/12.;
            break;
        
        case 4:
            a[0] = frac(lvl,u001,u101);
            b[0] = frac(lvl,u001,u011);
            c[0] = frac(lvl,u001,u000);
            a[1] = frac(lvl,u010,u110);
            b[1] = frac(lvl,u010,u000);
            c[1] = frac(lvl,u010,u011);
            deltav = (a[0]*b[0]*c[0]+a[1]*b[1]*c[1])/6.;
            break;
        
        case 5:
            a[0] = frac(lvl,u001,u101);
            a[1] = frac(lvl,u000,u100);
            a[2] = frac(lvl,u010,u110);
            b[0] = frac(lvl,u001,u011);
            c[0] = frac(lvl,u010,u011);
            deltav = (b[0]*a[0]+b[0]*a[1]+a[1]*(1.-c[0]*b[0])
                      +c[0]*a[1]+c[0]*a[2])/6.;
            break;
        
        case 6:
            a[0] = frac(lvl,u000,u100);
            a[1] = frac(lvl,u010,u110);
            a[2] = frac(lvl,u011,u111);
            a[3] = frac(lvl,u001,u101);
            deltav = (a[0]+a[1]+a[2]+a[3])/4.;
            break;
        
        case 7:
            a[0] = frac(lvl,u001,u101);
            b[0] = frac(lvl,u001,u011);
            c[0] = frac(lvl,u001,u000);
            a[1] = frac(lvl,u100,u000);
            b[1] = frac(lvl,u100,u110);
            c[1] = frac(lvl,u100,u101);
            a[2] = frac(lvl,u010,u110);
            b[2] = frac(lvl,u010,u000);
            c[2] = frac(lvl,u010,u011);
            deltav = (a[0]*b[0]*c[0]+a[1]*b[1]*c[1]+a[2]*b[2]*c[2])/6.;
            break;
        
        case 8:
            a[0] = frac(lvl,u001,u101);
            a[1] = frac(lvl,u010,u110);
            b[0] = frac(lvl,u001,u011);
            b[1] = frac(lvl,u100,u110);
            c[0] = frac(lvl,u010,u011);
            c[1] = frac(lvl,u100,u101);
            deltav = (a[0]*(b[0]+1.-b[1]*c[1])+a[1]*(c[0]*(1.-b[0])+1.+b[0]-b[1])
                      +b[1]*(1.+c[1])+1.)/6.;
            break;
        
        case 9:
            a[0] = frac(lvl,u011,u111);
            b[0] = frac(lvl,u011,u001);
            c[0] = frac(lvl,u011,u010);
            a[1] = frac(lvl,u100,u000);
            b[1] = frac(lvl,u100,u110);
            c[1] = frac(lvl,u100,u101);
            deltav = (a[0]*b[0]*c[0]+a[1]*b[1]*c[1])/6.;
            break;
        
        case 10:
            a[0] = frac(lvl,u011,u111);
            b[0] = frac(lvl,u011,u001);
            c[0] = frac(lvl,u011,u010);
            b[1] = frac(lvl,u000,u010);
            b[2] = frac(lvl,u100,u110);
            c[1] = frac(lvl,u000,u001);
            c[2] = frac(lvl,u100,u101);
            deltav = (2.*a[0]*b[0]*c[0]+b[1]*c[2]+2.*b[1]*c[1]+2.*b[2]*c[2]
                      +b[2]*c[1])/12.;
            break;
        
        case 11:
            a[0] = frac(lvl,u011,u111);
            a[1] = frac(lvl,u001,u101);
            b[0] = frac(lvl,u100,u110);
            b[1] = frac(lvl,u000,u010);
            c[0] = frac(lvl,u011,u010);
            c[1] = frac(lvl,u100,u101);
            deltav = (c[1]*b[0]+c[1]*b[1]+b[1]+b[1]*a[1]-a[1]*b[1]*c[1]
                      +a[1]-a[1]*b[1]*c[0]+a[0]*c[0]+a[1]*c[0])/6.;
            break;
        
        case 12:
            a[0] = frac(lvl,u011,u111);
            a[1] = frac(lvl,u010,u110);
            b[0] = frac(lvl,u011,u001);
            b[1] = frac(lvl,u100,u110);
            c[0] = frac(lvl,u000,u001);
            c[1] = frac(lvl,u100,u101);
            deltav = (a[0]*b[0]+a[1]*b[0]+a[1]+a[1]*c[0]-a[1]*b[0]*c[0]+c[0]
                      -a[1]*b[1]*c[0]+b[1]*c[1]+b[1]*c[0])/6.;
            break;
        
        case 13:
            a[0] = frac(lvl,u100,u000);
            a[1] = frac(lvl,u001,u101);
            a[2] = frac(lvl,u011,u111);
            a[3] = frac(lvl,u010,u110);
            b[0] = frac(lvl,u100,u110);
            b[1] = frac(lvl,u010,u000);
            c[0] = frac(lvl,u100,u101);
            c[1] = frac(lvl,u001,u000);
            deltav = (a[0]*b[0]*c[0]+c[1]*a[1]+c[1]*a[2]+a[2]*(1.-b[1]*c[1])
                      +b[1]*a[2]+b[1]*a[3])/6.;
            break;
        
        case 14:
            a[0] = frac(lvl,u011,u111);
            a[1] = frac(lvl,u010,u110);
            a[2] = frac(lvl,u101,u001);
            a[3] = frac(lvl,u100,u000);
            b[0] = frac(lvl,u011,u001);
            b[1] = frac(lvl,u010,u000);
            b[2] = frac(lvl,u101,u111);
            b[3] = frac(lvl,u100,u110);
            deltav = (a[0]*b[1]+2.*(a[0]*b[0]+a[1]*b[1])+a[1]*b[0]+a[2]*b[3]
                      +2.*(a[2]*b[2]+a[3]*b[3])+a[3]*b[2])/12.;
            break;
        
        case 15:
            a[0] = frac(lvl,u000,u100);
            b[0] = frac(lvl,u000,u010);
            c[0] = frac(lvl,u000,u001);
            a[1] = frac(lvl,u011,u111);
            b[1] = frac(lvl,u011,u001);
            c[1] = frac(lvl,u011,u010);
            a[2] = frac(lvl,u101,u001);
            b[2] = frac(lvl,u101,u111);
            c[2] = frac(lvl,u101,u100);
            a[3] = frac(lvl,u110,u010);
            b[3] = frac(lvl,u110,u100);
            c[3] = frac(lvl,u110,u111);
            deltav = (a[0]*b[0]*c[0]+a[1]*b[1]*c[1]
                      +a[2]*b[2]*c[2]+a[3]*b[3]*c[3])/6.;
            break;
        };
        return(deltav);
    }

    double CalcDv_(const double u000, const double u100, const double u110,
                   const double u010, const double u001, const double u101,
                   const double u111, const double u011, const double lvl,
                   const int sgn)
    {
        int             i, j, k;
        double  f[2][2][2];
        int             corn[2][2][2][3];
        int             val, temps;
        int             rot, newval, tempr[3];
        double  deltav;
        int             oval;
        static int lookup[128][2]
            = {0, 1, 0, 2, -2, 1, 0, 3, 3, 1, -1, 3, 0, 4, 0, 5, 1, 2, 1, 6,
               1, 3, 1, 7, 1, 10, -1, 7, 1, 11, 0, 6, -3, 1, 2, 3, -3, 9,
               -3, 13, 2, 9, 2, 11, 0, 7, 0, 8, 0, 9, 0, 10, -2, 37, 0, 11,
               -1, 37, 0, 12, 0, 13, -1, 47, -3, 2, -3, 6, -2, 3, -3, 14,
               1, 24, -3, 70, 1, 25, 1, 29, -2, 6, 1, 22, -2, 7, 1, 23, -2, 70,
               1, 30, 1, 27, -3, 206, -3, 3, -3, 7, -3, 11, -3, 15, 2, 25,
               2, 27, 3, 225, 3, 241, -2, 38, 3, 210, 3, 226, 3, 242, 0, 14,
               -3, 199, -3, 203, 3, 243, 2, 4, 2, 6, -1, 24, 2, 38, 3, 5,
               2, 14, -1, 25, -1, 27, 3, 6, -1, 22, 3, 38, -1, 30, 3, 7,
               3, 23, -1, 29, 3, 55, 2, 5, 2, 7, -3, 25, 2, 39, 2, 13, 2, 15,
               -2, 225, -2, 241, 3, 70, -1, 54, -1, 60, 1, 199, -2, 228,
               -2, 244, 1, 203, -2, 245, -3, 18, 2, 22, 2, 52, -3, 30, 2, 28,
               2, 30, 2, 60, -2, 211, 2, 148, 0, 15, 2, 180, 2, 182, -3, 210,
               -3, 214, -2, 199, -3, 222, 2, 21, 2, 23, -3, 27, 3, 244, 2, 29,
               2, 31, 2, 61, 3, 245, 1, 210, 1, 214, 3, 230, 3, 246, -2, 230,
               -2, 246, 3, 231, 3, 247
        };

        f[0][0][0] = u000;
        f[0][0][1] = u001;
        f[0][1][0] = u010;
        f[0][1][1] = u011;
        f[1][0][0] = u100;
        f[1][0][1] = u101;
        f[1][1][0] = u110;
        f[1][1][1] = u111;
        
        for (i=0; i<2; ++i)
            for (j=0; j<2; ++j)
                for (k=0; k<2; ++k) {
                    corn[i][j][k][0] = i;
                    corn[i][j][k][1] = j;
                    corn[i][j][k][2] = k;
                };
        val = sgn*(f[1][1][1]-lvl) > 0 ? 1 : 0;
        val = 2*val + (sgn*(f[1][1][0]-lvl) > 0 ? 1 : 0);
        val = 2*val + (sgn*(f[1][0][1]-lvl) > 0 ? 1 : 0);
        val = 2*val + (sgn*(f[1][0][0]-lvl) > 0 ? 1 : 0);
        val = 2*val + (sgn*(f[0][1][1]-lvl) > 0 ? 1 : 0);
        val = 2*val + (sgn*(f[0][1][0]-lvl) > 0 ? 1 : 0);
        val = 2*val + (sgn*(f[0][0][1]-lvl) > 0 ? 1 : 0);
        val = 2*val + (sgn*(f[0][0][0]-lvl) > 0 ? 1 : 0);
        oval = val;
        temps = 1;
        if (val > 127) {
            temps = -temps;
            val = 255-val;
        };
        rot = lookup[val][0];
        newval = lookup[val][1];
        
        while (rot != 0) {
            for (i=0; i<3; ++i) {
                switch(rot) {
                case 1:         /* x-plus rotation      */
                    tempr[i] = corn[0][0][0][i];
                    corn[0][0][0][i] = corn[0][0][1][i];
                    corn[0][0][1][i] = corn[0][1][1][i];
                    corn[0][1][1][i] = corn[0][1][0][i];
                    corn[0][1][0][i] = tempr[i];
                    tempr[i] = corn[1][0][0][i];
                    corn[1][0][0][i] = corn[1][0][1][i];
                    corn[1][0][1][i] = corn[1][1][1][i];
                    corn[1][1][1][i] = corn[1][1][0][i];
                    corn[1][1][0][i] = tempr[i];
                    break;
                case -1:        /* x-minus rotation     */
                    tempr[i] = corn[0][0][0][i];
                    corn[0][0][0][i] = corn[0][1][0][i];
                    corn[0][1][0][i] = corn[0][1][1][i];
                    corn[0][1][1][i] = corn[0][0][1][i];
                    corn[0][0][1][i] = tempr[i];
                    tempr[i] = corn[1][0][0][i];
                    corn[1][0][0][i] = corn[1][1][0][i];
                    corn[1][1][0][i] = corn[1][1][1][i];
                    corn[1][1][1][i] = corn[1][0][1][i];
                    corn[1][0][1][i] = tempr[i];
                    break;
                case 2:         /* y-plus rotation      */
                    tempr[i] = corn[0][0][0][i];
                    corn[0][0][0][i] = corn[1][0][0][i];
                    corn[1][0][0][i] = corn[1][0][1][i];
                    corn[1][0][1][i] = corn[0][0][1][i];
                    corn[0][0][1][i] = tempr[i];
                    tempr[i] = corn[0][1][0][i];
                    corn[0][1][0][i] = corn[1][1][0][i];
                    corn[1][1][0][i] = corn[1][1][1][i];
                    corn[1][1][1][i] = corn[0][1][1][i];
                    corn[0][1][1][i] = tempr[i];
                    break;
                case -2:        /* y-minus rotation     */
                    tempr[i] = corn[0][0][0][i];
                    corn[0][0][0][i] = corn[0][0][1][i];
                    corn[0][0][1][i] = corn[1][0][1][i];
                    corn[1][0][1][i] = corn[1][0][0][i];
                    corn[1][0][0][i] = tempr[i];
                    tempr[i] = corn[0][1][0][i];
                    corn[0][1][0][i] = corn[0][1][1][i];
                    corn[0][1][1][i] = corn[1][1][1][i];
                    corn[1][1][1][i] = corn[1][1][0][i];
                    corn[1][1][0][i] = tempr[i];
                    break;
                case 3:         /* z-plus rotation      */
                    tempr[i] = corn[0][0][0][i];
                    corn[0][0][0][i] = corn[0][1][0][i];
                    corn[0][1][0][i] = corn[1][1][0][i];
                    corn[1][1][0][i] = corn[1][0][0][i];
                    corn[1][0][0][i] = tempr[i];
                    tempr[i] = corn[0][0][1][i];
                    corn[0][0][1][i] = corn[0][1][1][i];
                    corn[0][1][1][i] = corn[1][1][1][i];
                    corn[1][1][1][i] = corn[1][0][1][i];
                    corn[1][0][1][i] = tempr[i];
                    break;
                case -3:        /* z-minus rotation     */
                    tempr[i] = corn[0][0][0][i];
                    corn[0][0][0][i] = corn[1][0][0][i];
                    corn[1][0][0][i] = corn[1][1][0][i];
                    corn[1][1][0][i] = corn[0][1][0][i];
                    corn[0][1][0][i] = tempr[i];
                    tempr[i] = corn[0][0][1][i];
                    corn[0][0][1][i] = corn[1][0][1][i];
                    corn[1][0][1][i] = corn[1][1][1][i];
                    corn[1][1][1][i] = corn[0][1][1][i];
                    corn[0][1][1][i] = tempr[i];
                    break;
                };
            };
            val = newval;
            if (val > 127) {
                temps = -temps;
                val = 255-val;
            };
            rot = lookup[val][0];
            newval = lookup[val][1];
                
        };
        val = newval;
        deltav = CalcRot_(f[corn[0][0][0][0]][corn[0][0][0][1]][corn[0][0][0][2]],
                          f[corn[1][0][0][0]][corn[1][0][0][1]][corn[1][0][0][2]],
                          f[corn[1][1][0][0]][corn[1][1][0][1]][corn[1][1][0][2]],
                          f[corn[0][1][0][0]][corn[0][1][0][1]][corn[0][1][0][2]],
                          f[corn[0][0][1][0]][corn[0][0][1][1]][corn[0][0][1][2]],
                          f[corn[1][0][1][0]][corn[1][0][1][1]][corn[1][0][1][2]],
                          f[corn[1][1][1][0]][corn[1][1][1][1]][corn[1][1][1][2]],
                          f[corn[0][1][1][0]][corn[0][1][1][1]][corn[0][1][1][2]],
                          lvl,val);
        
        if (temps < 0)
            deltav = 1.-deltav;
        return(deltav);
    }

    double UniformMesh3D::Volume(const int l, const double val) 
    {
        double vol = 0.;
   
        for (int i=0; i<maxi-1; ++i)
            for (int j=0; j<maxj-1; ++j)
                for (int k=0; k<maxk-1; ++k)
                    vol += dx*dy*dz*CalcDv_(data_(i,j,k,l), data_(i+1,j,k,l),
                                            data_(i+1,j+1,k,l), data_(i,j+1,k,l),
                                            data_(i,j,k+1,l), data_(i+1,j,k+1,l),
                                            data_(i+1,j+1,k+1,l), data_(i,j+1,k+1,l),
                                            val, 1);
        return(vol);
    }

    char UniformMesh3D::IsCrossing(const double value, const double y0,
                                   const double y1, double& f)
    {
        f = frac(value, y0, y1);
        return f>=0 && f<1;
    }

    char UniformMesh3D::TestValidity(const char* file, const int line,
                                     const int kind) const
    {
        char ans = true;
   
        if (kind >= 0) {
            for (int iii=0; iii<maxi; ++iii)
                for (int jjj=0; jjj<maxj; ++jjj) 
                    for (int kkk=0; kkk<maxk; ++kkk) {
                        if (!finite(data_(iii,jjj,kkk,kind))) {
                            std::cerr << "Warning: Invalid data at (" << iii << ',' << jjj
                                      << ',' << kkk << ',' << "kFunc) in " << file
                                      << ':' << line << '\n';
                            ans = false;
                        }
                    }
        } else {
            for (int l=0; l<maxl; ++l) 
                for (int iii=0; iii<maxi; ++iii)
                    for (int jjj=0; jjj<maxj; ++jjj) 
                        for (int kkk=0; kkk<maxk; ++kkk) {
                            if (!finite(data_(iii,jjj,kkk,l))) {
                                std::cerr << "Warning: Invalid data at (" << iii << ',' << jjj
                                          << ',' << kkk << ',' << l << ") in " << file
                                          << ':' << line << '\n';
                                ans = false;
                            }
                        }
        }

        if (ans) std::cerr << "Data OK in " << file << ':' << line << '\n';
        return ans;   
    }

    void UniformMesh3D::MinMax(const int l, double& themin, double& themax) const
    {
        themin = data_(0,0,0,l);
        themax = data_(0,0,0,l);
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j) 
                for (int k=0; k<maxk; ++k) {
                    if (data_(i,j,k,l) < themin) themin = data_(i,j,k,l);
                    if (data_(i,j,k,l) > themax) themax = data_(i,j,k,l); 
                }
    }
        
}

