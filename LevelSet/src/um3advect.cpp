/****************************************************************
    $Header: um3advect.cpp,v 2.3 99/01/06 13:59:06 chopp Exp $
  
    $Log:       um3advect.cpp,v $
    Revision 2.3  99/01/06  13:59:06  13:59:06  chopp (David Chopp)
    *** none ***

    Revision 1.7  98/05/18  14:06:12  14:06:12  chopp (David Chopp)
    *** empty log message ***

    ***********************************************************/
#include "numtrait.h"
#include "datavec.h"
#include "uniformmesh3d.h"
#include "um3heap.h"
#include "debug.h"
#include "plotwindow3d.h"
#include "um3boundary.h"
#include "utility.h"
#include "tricubic.h"
#ifdef MEMWATCH
#include "memwatch.h"
#endif

namespace levelset {
        
    void UniformMesh3D::Advance(const int l, const int kv, const int knorm,
                                const double dt)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,l) += dt*data_(i,j,k,knorm)*data_(i,j,k,kv);
        bc->Apply(l);
    }

    void UniformMesh3D::Advance(const int l, const int kv, const double dt)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,l) += dt*data_(i,j,k,kv);
        bc->Apply(l);
    }

    void UniformMesh3D::Product(const int l, const int la, const int lb)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,l) = data_(i,j,k,la)*data_(i,j,k,lb);
        bc->Apply(k);
    }

    double UniformMesh3D::ComputeTVal(const int i, const int j, const int k, 
                                      const int l)
    {
        double a, b, c, d;
   
        if (data_(i,j,k,l) > 0) {
            if (i <= 0)
                a = data_(i+1,j,k,l);
            else if (i >= maxi-1)
                a = data_(i-1,j,k,l);
            else
                a = min(data_(i-1,j,k,l),data_(i+1,j,k,l));

            if (j <= 0)
                b = data_(i,j+1,k,l);
            else if (j >= maxj-1)
                b = data_(i,j-1,k,l);
            else
                b = min(data_(i,j-1,k,l),data_(i,j+1,k,l));

            if (k <= 0)
                c = data_(i,j,k+1,l);
            else if (k >= maxk-1)
                c = data_(i,j,k-1,l);
            else
                c = min(data_(i,j,k-1,l),data_(i,j,k+1,l));

            data_(i,j,k,l) = min(a+dx,b+dy,c+dz);
            if ((d=dx*dx+dy*dy-sqr(a-b)) >= 0.) 
                data_(i,j,k,l) = min(data_(i,j,k,l),
                                     (a*dy*dy+b*dx*dx+dx*dy*sqrt(d))/(dx*dx+dy*dy));
            if ((d=dx*dx+dz*dz-sqr(a-c)) >= 0.) 
                data_(i,j,k,l) = min(data_(i,j,k,l),
                                     (a*dz*dz+c*dx*dx+dx*dz*sqrt(d))/(dx*dx+dz*dz));
            if ((d=dy*dy+dz*dz-sqr(b-c)) >= 0.) 
                data_(i,j,k,l) = min(data_(i,j,k,l),
                                     (b*dz*dz+c*dy*dy+dy*dz*sqrt(d))/(dy*dy+dz*dz));
            if ((d=1./dx/dx+1./dy/dy+1./dz/dz-sqr((a-b)/dx/dy)-sqr((a-c)/dx/dz)
                 -sqr((b-c)/dy/dz)) >= 0.)
                data_(i,j,k,l) = min(data_(i,j,k,l),
                                     (a/dx/dx+b/dy/dy+c/dz/dz+sqrt(d)
                                         )/(1./dx/dx+1./dy/dy+1./dz/dz));
        }
        else {
            if (i <= 0)
                a = data_(i+1,j,k,l);
            else if (i >= maxi-1)
                a = data_(i-1,j,k,l);
            else
                a = max(data_(i-1,j,k,l),data_(i+1,j,k,l));

            if (j <= 0)
                b = data_(i,j+1,k,l);
            else if (j >= maxj-1)
                b = data_(i,j-1,k,l);
            else
                b = max(data_(i,j-1,k,l),data_(i,j+1,k,l));

            if (k <= 0)
                c = data_(i,j,k+1,l);
            else if (k >= maxk-1)
                c = data_(i,j,k-1,l);
            else
                c = max(data_(i,j,k-1,l),data_(i,j,k+1,l));

            data_(i,j,k,l) = max(a-dx,b-dy,c-dz);
            if ((d=dx*dx+dy*dy-sqr(a-b)) >= 0.) 
                data_(i,j,k,l) = max(data_(i,j,k,l),
                                     (a*dy*dy+b*dx*dx-dx*dy*sqrt(d))/(dx*dx+dy*dy));
            if ((d=dx*dx+dz*dz-sqr(a-c)) >= 0.) 
                data_(i,j,k,l) = max(data_(i,j,k,l),
                                     (a*dz*dz+c*dx*dx-dx*dz*sqrt(d))/(dx*dx+dz*dz));
            if ((d=dy*dy+dz*dz-sqr(b-c)) >= 0.) 
                data_(i,j,k,l) = max(data_(i,j,k,l),
                                     (b*dz*dz+c*dy*dy-dy*dz*sqrt(d))/(dy*dy+dz*dz));
            if ((d=1./dx/dx+1./dy/dy+1./dz/dz-sqr((a-b)/dx/dy)-sqr((a-c)/dx/dz)
                 -sqr((b-c)/dy/dz)) >= 0.)
                data_(i,j,k,l) = max(data_(i,j,k,l),
                                     (a/dx/dx+b/dy/dy+c/dz/dz-sqrt(d)
                                         )/(1./dx/dx+1./dy/dy+1./dz/dz));
        }
      
        return data_(i,j,k,l);
    }

    void UniformMesh3D::ExtendVals(const int i, const int j, const int k,
                                   const int l, const int ikus, 
                                   const int* p, const int pl)
    {
        int d;
        int sgn = data_(i,j,k,l) >= 0. ? 1 : -1;
        d  = ((i > 0) && idata_(i-1,j,k,ikus)==2) ? 1 : 0;
        d += ((j > 0) && idata_(i,j-1,k,ikus)==2) ? 2 : 0;
        d += ((k > 0) && idata_(i,j,k-1,ikus)==2) ? 4 : 0;
        d += ((i < maxi-1) && idata_(i+1,j,k,ikus)==2) ? 8 : 0;
        d += ((j < maxj-1) && idata_(i,j+1,k,ikus)==2) ? 16 : 0;
        d += ((k < maxk-1) && idata_(i,j,k+1,ikus)==2) ? 32 : 0;

        if ((d & 9) == 9) {
            if (sgn*data_(i-1,j,k,l) < sgn*data_(i+1,j,k,l))
                d -= 8;
            else
                d -= 1;
        }
        if ((d & 18) == 18) {
            if (sgn*data_(i,j-1,k,l) < sgn*data_(i,j+1,k,l))
                d -= 16;
            else
                d -= 2;
        }
        if ((d & 36) == 36) {
            if (sgn*data_(i,j,k-1,l) < sgn*data_(i,j,k+1,l))
                d -= 32;
            else
                d -= 4;
        }
      
        int m;
        double a = 0.;
        double b = 0.;
        switch(d) {
        case 0:
            std::cout << "Case 0 hit\n";
            for (m=0; m<pl; ++m) {
                if (idata_(i+1,j,k,ikus)) {
                    a += data_(i+1,j,k,p[m]);
                    b += 1.;
                }
                if (idata_(i-1,j,k,ikus)) {
                    a += data_(i-1,j,k,p[m]);
                    b += 1.;
                }
                if (idata_(i,j+1,k,ikus)) {
                    a += data_(i,j+1,k,p[m]);
                    b += 1.;
                }
                if (idata_(i,j-1,k,ikus)) {
                    a += data_(i,j-1,k,p[m]);
                    b += 1.;
                }
                if (idata_(i,j,k-1,ikus)) {
                    a += data_(i,j,k-1,p[m]);
                    b += 1.;
                }
                if (idata_(i,j,k+1,ikus)) {
                    a += data_(i,j,k+1,p[m]);
                    b += 1.;
                }
                data_(i,j,k,p[m]) = a/b;
            }
            break;
        case 1:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = data_(i-1,j,k,p[m]);
            break;
        case 2:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = data_(i,j-1,k,p[m]);
            break;
        case 3:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j-1,k,l))*data_(i,j-1,k,p[m])
                                     +dy*dy*(data_(i,j,k,l)-data_(i-1,j,k,l))
                                     *data_(i-1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j-1,k,l))
                      +dy*dy*(data_(i,j,k,l)-data_(i-1,j,k,l)));
            break;
        case 4:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = data_(i,j,k-1,p[m]);
            break;
        case 5:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j,k-1,l))*data_(i,j,k-1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i-1,j,k,l))
                                     *data_(i-1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j,k-1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i-1,j,k,l)));
            break;      
        case 6:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dy*dy*(data_(i,j,k,l)
                                            -data_(i,j,k-1,l))*data_(i,j,k-1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i,j-1,k,l))
                                     *data_(i,j-1,k,p[m]))
                    /(dy*dy*(data_(i,j,k,l)-data_(i,j,k-1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i,j-1,k,l)));
            break;
        case 7:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i-1,j,k,l)
                                         )*data_(i-1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j-1,k,l)
                                         )*data_(i,j-1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k-1,l)
                                         )*data_(i,j,k-1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i-1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j-1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k-1,l))/dz/dz);
            break;
        case 8:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = data_(i+1,j,k,p[m]);
            break;
        case 10:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j-1,k,l))*data_(i,j-1,k,p[m])
                                     +dy*dy*(data_(i,j,k,l)-data_(i+1,j,k,l))
                                     *data_(i+1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j-1,k,l))
                      +dy*dy*(data_(i,j,k,l)-data_(i+1,j,k,l)));
            break;
        case 12:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j,k-1,l))*data_(i,j,k-1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i+1,j,k,l))
                                     *data_(i+1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j,k-1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i+1,j,k,l)));
            break;      
        case 14:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i+1,j,k,l)
                                         )*data_(i+1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j-1,k,l)
                                         )*data_(i,j-1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k-1,l)
                                         )*data_(i,j,k-1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i+1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j-1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k-1,l))/dz/dz);
            break;
        case 16:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = data_(i,j+1,k,p[m]);
            break;
        case 17:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j+1,k,l))*data_(i,j+1,k,p[m])
                                     +dy*dy*(data_(i,j,k,l)-data_(i-1,j,k,l))
                                     *data_(i-1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j+1,k,l))
                      +dy*dy*(data_(i,j,k,l)-data_(i-1,j,k,l)));
            break;
        case 20:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dy*dy*(data_(i,j,k,l)
                                            -data_(i,j,k-1,l))*data_(i,j,k-1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i,j+1,k,l))
                                     *data_(i,j+1,k,p[m]))
                    /(dy*dy*(data_(i,j,k,l)-data_(i,j,k-1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i,j+1,k,l)));
            break;
        case 21:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i-1,j,k,l)
                                         )*data_(i-1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j+1,k,l)
                                         )*data_(i,j+1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k-1,l)
                                         )*data_(i,j,k-1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i-1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j+1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k-1,l))/dz/dz);
            break;
        case 24:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j+1,k,l))*data_(i,j+1,k,p[m])
                                     +dy*dy*(data_(i,j,k,l)-data_(i+1,j,k,l))
                                     *data_(i+1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j+1,k,l))
                      +dy*dy*(data_(i,j,k,l)-data_(i+1,j,k,l)));
            break;
        case 28:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i+1,j,k,l)
                                         )*data_(i+1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j+1,k,l)
                                         )*data_(i,j+1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k-1,l)
                                         )*data_(i,j,k-1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i+1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j+1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k-1,l))/dz/dz);
            break;
        case 32:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = data_(i,j,k+1,p[m]);
            break;
        case 33:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j,k+1,l))*data_(i,j,k+1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i-1,j,k,l))
                                     *data_(i-1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j,k+1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i-1,j,k,l)));
            break;      
        case 34:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dy*dy*(data_(i,j,k,l)
                                            -data_(i,j,k+1,l))*data_(i,j,k+1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i,j-1,k,l))
                                     *data_(i,j-1,k,p[m]))
                    /(dy*dy*(data_(i,j,k,l)-data_(i,j,k+1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i,j-1,k,l)));
            break;
        case 35:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i-1,j,k,l)
                                         )*data_(i-1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j-1,k,l)
                                         )*data_(i,j-1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k+1,l)
                                         )*data_(i,j,k+1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i-1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j-1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k+1,l))/dz/dz);
            break;
        case 40:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dx*dx*(data_(i,j,k,l)
                                            -data_(i,j,k+1,l))*data_(i,j,k+1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i+1,j,k,l))
                                     *data_(i+1,j,k,p[m]))
                    /(dx*dx*(data_(i,j,k,l)-data_(i,j,k+1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i+1,j,k,l)));
            break;
        case 42:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i+1,j,k,l)
                                         )*data_(i+1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j-1,k,l)
                                         )*data_(i,j-1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k+1,l)
                                         )*data_(i,j,k+1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i+1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j-1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k+1,l))/dz/dz);
            break;
        case 48:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = (dy*dy*(data_(i,j,k,l)
                                            -data_(i,j,k+1,l))*data_(i,j,k+1,p[m])
                                     +dz*dz*(data_(i,j,k,l)-data_(i,j+1,k,l))
                                     *data_(i,j+1,k,p[m]))
                    /(dy*dy*(data_(i,j,k,l)-data_(i,j,k+1,l))
                      +dz*dz*(data_(i,j,k,l)-data_(i,j+1,k,l)));
            break;
        case 49:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i-1,j,k,l)
                                         )*data_(i-1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j+1,k,l)
                                         )*data_(i,j+1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k+1,l)
                                         )*data_(i,j,k+1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i-1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j+1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k+1,l))/dz/dz);
            break;
        case 56:
            for (m=0; m<pl; ++m)
                data_(i,j,k,p[m]) = ((data_(i,j,k,l)-data_(i+1,j,k,l)
                                         )*data_(i+1,j,k,p[m])/dx/dx
                                     +(data_(i,j,k,l)-data_(i,j+1,k,l)
                                         )*data_(i,j+1,k,p[m])/dy/dy
                                     +(data_(i,j,k,l)-data_(i,j,k+1,l)
                                         )*data_(i,j,k+1,p[m])/dz/dz)
                    /((data_(i,j,k,l)-data_(i+1,j,k,l))/dx/dx
                      +(data_(i,j,k,l)-data_(i,j+1,k,l))/dy/dy
                      +(data_(i,j,k,l)-data_(i,j,k+1,l))/dz/dz);
            break;
        default:
            std::cerr << "Oops! in ExtendVals, case = " << d << '\n';
        }
    }

    void UniformMesh3D::ExtendFunc(const int i, const int j, const int k,
                                   const int l, const int ikus)
    {
        int sgn = data_(i,j,k,l) >= 0. ? 1 : -1;

        double d = 0;
        data_(i,j,k,l) = 0.;
   
        if (i>0 && idata_(i-1,j,k,ikus) == 2) {
            if (i>1 && idata_(i-2,j,k,ikus) == 2)
                data_(i,j,k,l) += (2*data_(i-1,j,k,l)-data_(i-2,j,k,l))/dx/dx;
            else
                data_(i,j,k,l) += (sgn*dx+data_(i-1,j,k,l))/dx/dx;
            d += 1./dx/dx;
        }
        if (i<maxi-1 && idata_(i+1,j,k,ikus) == 2) {
            if (i<maxi-2 && idata_(i+2,j,k,ikus) == 2)
                data_(i,j,k,l) += (2*data_(i+1,j,k,l)-data_(i+2,j,k,l))/dx/dx;
            else
                data_(i,j,k,l) += (sgn*dx+data_(i+1,j,k,l))/dx/dx;
            d += 1./dx/dx;
        }
        if (j>0 && idata_(i,j-1,k,ikus) == 2) {
            if (j>1 && idata_(i,j-2,k,ikus) == 2)
                data_(i,j,k,l) += 2*data_(i,j-1,k,l)-data_(i,j-2,k,l);
            else
                data_(i,j,k,l) += sgn*dy+data_(i,j-1,k,l);
            d += 1./dy/dy;
        }
        if (j<maxj-1 && idata_(i,j+1,k,ikus) == 2) {
            if (j<maxj-2 && idata_(i,j+2,k,ikus) == 2)
                data_(i,j,k,l) += (2*data_(i,j+1,k,l)-data_(i,j+2,k,l))/dy/dy;
            else
                data_(i,j,k,l) += (sgn*dy+data_(i,j+1,k,l))/dy/dy;
            d += 1./dy/dy;
        }
        if (k>0 && idata_(i,j,k-1,ikus) == 2) {
            if (k>1 && idata_(i,j,k-2,ikus) == 2)
                data_(i,j,k,l) += 2*data_(i,j,k-1,l)-data_(i,j,k-2,l);
            else
                data_(i,j,k,l) += sgn*dz+data_(i,j,k-1,l);
            d += 1./dz/dz;
        }
        if (k<maxk-1 && idata_(i,j,k+1,ikus) == 2) {
            if (k<maxk-2 && idata_(i,j,k+2,ikus) == 2)
                data_(i,j,k,l) += (2*data_(i,j,k+1,l)-data_(i,j,k+2,l))/dz/dz;
            else
                data_(i,j,k,l) += (sgn*dz+data_(i,j,k+1,l))/dz/dz;
            d += 1./dz/dz;
        }

        data_(i,j,k,l) /= d;
    }

#if 0
    void UniformMesh3D::UpdateNeighbor(UniformMesh3D_Heap& heap, const int i,
                                       const int j, const int k, const int l, 
                                       const int ikus, const int ikhi)
    {
        switch(idata_(i,j,k,ikus)) {
        case 0:
            heap.Insert(i,j,k,abs(ComputeTVal(i,j,k,l)),ikhi);
            idata_(i,j,k,ikus) = 1;
            break;
        case 1:
            heap.Change(idata_(i,j,k,ikhi),abs(ComputeTVal(i,j,k,l)),ikhi);
            break;
        case 2:
            break;
        default:
            cerr << "Error: Gridpoint " << i << "," << j << "," << k 
                 << " not initialized.\n";
            break;
        }
    }
#else
    void UniformMesh3D::UpdateNeighbor(UniformMesh3D_Heap& heap, const int i,
                                       const int j, const int k, const int l, 
                                       const int ikus, const int ikhi)
    {
        switch(idata_(i,j,k,ikus)) {
        case 0:
            heap.Insert(i,j,k,fabs(ComputeTVal(i,j,k,l,ikus)),ikhi);
            idata_(i,j,k,ikus) = 1;
            break;
        case 1:
            heap.Change(idata_(i,j,k,ikhi),fabs(ComputeTVal(i,j,k,l,ikus)),ikhi);
            break;
        case 2:
            break;
        default:
            std::cerr << "Error: Gridpoint " << i << "," << j << "," << k 
                      << " not initialized.\n";
            break;
        }
    }
#endif

#if 1
    void UniformMesh3D::InitLevelSet(const int i, const int j, const int k, const int func, 
                                     const int phi, const int ikus)
    {
        Tricubic p;
        double temp[4][4][4];
        for (int ii=0; ii<4; ++ii)
            for (int jj=0; jj<4; ++jj)
                for (int kk=0; kk<4; ++kk)
                    temp[ii][jj][kk] = data_(i-1+ii,j-1+jj,k-1+kk,func);
        p.BuildwDeriv(temp, dx, dy, dz);
        
        double ax, ay, az;
        double dist[2][2][2];
        char clean[2][2][2];
        for (int ii=0; ii<2; ++ii)
            for (int jj=0; jj<2; ++jj)
                for (int kk=0; kk<2; ++kk) 
                    dist[ii][jj][kk] = p.LocalDist(ii*dx, jj*dy, kk*dz, ax, ay, az, clean[ii][jj][kk]);
   
        for (int ii=0; ii<2; ++ii)
            for (int jj=0; jj<2; ++jj) 
                for (int kk=0; kk<2; ++kk) {
                    if (clean[ii][jj][kk]) {
                        if (idata_(i+ii,j+jj,k+kk,ikus) < 2) {
                            data_(i+ii,j+jj,k+kk,phi) = copysign(dist[ii][jj][kk],
                                                                 data_(i+ii,j+jj,k+kk,func));
                            idata_(i+ii,j+jj,k+kk,ikus) = 2;
                        }
                        else if (fabs(data_(i+ii,j+jj,k+kk,phi)) > dist[ii][jj][kk]) {
                            data_(i+ii,j+jj,k+kk,phi) =copysign(dist[ii][jj][kk],
                                                                data_(i+ii,j+jj,k+kk,func));
                            idata_(i+ii,j+jj,k+kk,ikus) = 2;
                        }
                    } else if (fabs(data_(i+ii,j+jj,k+kk,phi)) > dist[ii][jj][kk]) {
                        data_(i+ii,j+jj,k+kk,phi) =copysign(dist[ii][jj][kk],data_(i+ii,j+jj,k+kk,func));
                        idata_(i+ii,j+jj,k+kk,ikus) = 1;
                    }
                }
    }   

#else
    void UniformMesh3D::InitLevelSet(const int i, const int j, const int k,
                                     const int func, const int phi, const int ikus,
                                     const int ikhi)
    {
        if (idata_(i,j,k,ikhi) == idata_(i,j,k,ikus)) {

            double s[3][2];

            if (idata_(i,j,k,ikus) & DLeft) 
                s[0][0] = frac(0., data_(i,j,k,func), data_(i-1,j,k,func));
            if (idata_(i,j,k,ikus) & DRight) 
                s[0][1] = frac(0., data_(i,j,k,func), data_(i+1,j,k,func));
            if (idata_(i,j,k,ikus) & DDown) 
                s[1][0] = frac(0., data_(i,j,k,func), data_(i,j-1,k,func));
            if (idata_(i,j,k,ikus) & DUp)
                s[1][1] = frac(0., data_(i,j,k,func), data_(i,j+1,k,func));
            if (idata_(i,j,k,ikus) & DBack) 
                s[2][0] = frac(0., data_(i,j,k,func), data_(i,j,k-1,func));
            if (idata_(i,j,k,ikus) & DAhead)
                s[2][1] = frac(0., data_(i,j,k,func), data_(i,j,k+1,func));      
            int sgn = (data_(i,j,k,func) >= 0) ? 1 : -1;
        
            int d = idata_(i,j,k,ikus);
            if ((d & DLeft) && (d & DRight)) 
                d -= s[0][0] < s[0][1] ? DRight : DLeft;
            if ((d & DDown) && (d & DUp))
                d -= s[1][0] < s[1][1] ? DUp : DDown;
            if ((d & DBack) && (d & DAhead))
                d -= s[2][0] < s[2][1] ? DAhead : DBack;
      
            switch(d) {
            case 0: break;
            case DLeft: data_(i,j,k,phi) = sgn*s[0][0]*dx; break;
            case DDown: data_(i,j,k,phi) = sgn*s[1][0]*dy; break;
            case DLeft | DDown:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][0]*dx)+1./sqr(s[1][0]*dy));
                break;
            case DRight: data_(i,j,k,phi) = sgn*s[0][1]*dx; break;
            case DRight | DDown:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[0][1]*dx));
                break;
            case DUp: data_(i,j,k,phi) = sgn*s[1][1]*dy; break;
            case DLeft | DUp:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][0]*dx)+1./sqr(s[1][1]*dy));
                break;
            case DUp | DRight:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][1]*dx)+1./sqr(s[1][1]*dy));
                break;
            case DBack:
                data_(i,j,k,phi) = sgn*s[2][0]*dz;
                break;
            case DLeft | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][0]*dx)+1./sqr(s[2][0]*dz));
                break;
            case DDown | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[2][0]*dz));
                break;
            case DLeft | DDown | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[2][0]*dz)
                                            +1./sqr(s[0][0]*dx));
                break;
            case DRight | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][1]*dx)+1./sqr(s[2][0]*dz));
                break;
            case DDown | DRight | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[2][0]*dz)
                                            +1./sqr(s[0][1]*dx));
                break;
            case DUp | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][1]*dy)+1./sqr(s[2][0]*dz));
                break;
            case DLeft | DUp | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][1]*dy)+1./sqr(s[2][0]*dz)
                                            +1./sqr(s[0][0]*dx));
                break;
            case DRight | DUp | DBack:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][1]*dy)+1./sqr(s[2][0]*dz)
                                            +1./sqr(s[0][1]*dx));
                break;
            case DAhead:
                data_(i,j,k,phi) = sgn*s[2][1]*dz;
                break;
            case DLeft | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][0]*dx)+1./sqr(s[2][1]*dz));
                break;
            case DDown | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[2][1]*dz));
                break;
            case DLeft | DDown | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[2][1]*dz)
                                            +1./sqr(s[0][0]*dx));
                break;
            case DRight | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[0][1]*dx)+1./sqr(s[2][1]*dz));
                break;
            case DDown | DRight | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][0]*dy)+1./sqr(s[2][1]*dz)
                                            +1./sqr(s[0][1]*dx));
                break;
            case DUp | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][1]*dy)+1./sqr(s[2][1]*dz));
                break;
            case DLeft | DUp | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][1]*dy)+1./sqr(s[2][1]*dz)
                                            +1./sqr(s[0][0]*dx));
                break;
            case DRight | DUp | DAhead:
                data_(i,j,k,phi) = sgn/sqrt(1./sqr(s[1][1]*dy)+1./sqr(s[2][1]*dz)
                                            +1./sqr(s[0][1]*dx));
                break;
            }
        }
    }
#endif

    void UniformMesh3D::ExtendVelocity(const int ktemp, const int kvel,
                                       const int ikus, const int ikhi,
                                       const int imask)
    {  
        NumericTrait<double> t;
      
        // Initialize the accepted and tentative points
      
        // Initialize the UPDATE_STATUS

        int i, j, k;

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    if (idata_(i,j,k,imask)) {
                        idata_(i,j,k,ikus) = 2;
                        data_(i,j,k,ktemp) = 0.;
                    } else {
                        idata_(i,j,k,ikus) = 0;
                        data_(i,j,k,ktemp) = t.max_value;
                    }
   
        // Now do fast marching to complete extension

        UniformMesh3D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    if (idata_(i,j,k,ikus) == 2) {
                        if (i > 0)
                            UpdateNeighbor(theHeap, i-1, j, k, ktemp, ikus, ikhi);
                        if (j > 0)
                            UpdateNeighbor(theHeap, i, j-1, k, ktemp, ikus, ikhi);
                        if (k > 0)
                            UpdateNeighbor(theHeap, i, j, k-1, ktemp, ikus, ikhi);
                        if (i < maxi-1)
                            UpdateNeighbor(theHeap, i+1, j, k, ktemp, ikus, ikhi);
                        if (j < maxj-1)
                            UpdateNeighbor(theHeap, i, j+1, k, ktemp, ikus, ikhi);
                        if (k < maxk-1)
                            UpdateNeighbor(theHeap, i, j, k+1, ktemp, ikus, ikhi);
                    }
            
        {
            int i, j, k;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,k,ikhi);
                idata_(i,j,k,ikus) = 2;
                if (!idata_(i,j,k,imask))
                    ExtendVals(i, j, k, ktemp, ikus, &kvel, 1);

                if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, ktemp, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, ktemp, ikus, ikhi);
                if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, ktemp, ikus, ikhi);
                if (i < maxi-1)
                    UpdateNeighbor(theHeap, i+1, j, k, ktemp, ikus, ikhi);
                if (j < maxj-1)
                    UpdateNeighbor(theHeap, i, j+1, k, ktemp, ikus, ikhi);
                if (k < maxk-1)
                    UpdateNeighbor(theHeap, i, j, k+1, ktemp, ikus, ikhi);
            }
        }

    }      

#if 1
    void UniformMesh3D::Reinitialize(const int func, const int phi, const int ikus,
                                     const int ikhi, const int dir)
    {
        NumericTrait<double> t;

        // Initialize the accepted and tentative points

        // Initialize the ikus

        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                for (k=0; k<maxk; ++k) {
                    idata_(i,j,k,ikus) = 0;
                    idata_(i,j,k,ikhi) = -1;
                    data_(i,j,k,phi) = data_(i,j,k,func) >= 0. ? t.max_value : t.min_value;
                }

        int cross;
        for (i=0; i<maxi-1; ++i)
            for (j=0; j<maxj-1; ++j) 
                for (k=0; k<maxk-1; ++k) {
                    cross = data_(i,j,k,func) > 0. ? 1 : 0;
                    cross += data_(i+1,j,k,func) > 0 ? 2 : 0;
                    cross += data_(i,j+1,k,func) > 0 ? 4 : 0;
                    cross += data_(i+1,j+1,k,func) > 0 ? 8 : 0;
                    cross += data_(i,j,k+1,func) > 0 ? 16 : 0;
                    cross += data_(i+1,j,k+1,func) > 0 ? 32 : 0;
                    cross += data_(i,j+1,k+1,func) > 0 ? 64 : 0;
                    cross += data_(i+1,j+1,k+1,func) > 0 ? 128 : 0;
                    if (cross > 0 && cross < 255) 
                        InitLevelSet(i,j,k,func,phi,ikus);
                }
                
        // Now do fast marching to complete extension

        UniformMesh3D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                for (k=0; k<maxk; ++k) {
                    if (idata_(i,j,k,ikus)==2) {
                        if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, phi, ikus, ikhi);
                        if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, phi, ikus, ikhi);
                        if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, phi, ikus, ikhi);
                        if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, k, phi, ikus, ikhi);
                        if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, k, phi, ikus, ikhi);
                        if (k < maxk-1) UpdateNeighbor(theHeap, i, j, k+1, phi, ikus, ikhi);
                    }
                }

        {
            int i, j;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,k,ikhi);
                idata_(i,j,k,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, phi, ikus, ikhi);
                if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, phi, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, k, phi, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, k, phi, ikus, ikhi);
                if (k < maxk-1) UpdateNeighbor(theHeap, i, j, k+1, phi, ikus, ikhi);
            }
        }

        switch (dir) {
        case 0:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    for (k=0; k<maxk; ++k)
                        data_(i,j,k,func) = data_(i,j,k,phi);
            break;
        case 1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    for (k=0; k<maxk; ++k)
                        if (data_(i,j,k,phi) > 0) data_(i,j,k,func) = data_(i,j,k,phi);
            break;
        case -1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    for (k=0; k<maxk; ++k)
                        if (data_(i,j,k,phi) < 0) data_(i,j,k,func) = data_(i,j,k,phi);
            break;
        }

        bc->Apply(func);

    }

#else
    void UniformMesh3D::Reinitialize(const int func, const int phi, const int ikus,
                                     const int ikhi, const int dir)
    {
        NumericTrait<double> t;

        // Initialize the accepted and tentative points

        // Initialize the ikus

        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    idata_(i,j,k,ikus) = 0;
                    idata_(i,j,k,ikhi) = 0;
                    data_(i,j,k,phi) = (data_(i,j,k,func) >= 0. ?
                                        t.max_value : t.min_value);
                }

        int cross;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {

                    if (i < maxi-1 && data_(i,j,k,func)*data_(i+1,j,k,func) <= 0) {
                        idata_(i,j,k,ikus) += DRight;
                        idata_(i+1,j,k,ikus) += DLeft;
                    }

                    if (j < maxj-1 && data_(i,j,k,func)*data_(i,j+1,k,func) <= 0) {
                        idata_(i,j,k,ikus) += DUp;
                        idata_(i,j+1,k,ikus) += DDown;
                    }

                    if (k < maxk-1 && data_(i,j,k,func)*data_(i,j,k+1,func) <= 0) {
                        idata_(i,j,k,ikus) += DAhead;
                        idata_(i,j,k+1,ikus) += DBack;
                    }
                }


        // Now do fast marching to complete extension

        UniformMesh3D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    idata_(i,j,k,ikhi) = idata_(i,j,k,ikus);
                    if (idata_(i,j,k,ikus)) {
                        InitLevelSet(i,j,k,func,phi,ikus,ikhi);
                        idata_(i,j,k,ikus) = 2;
                    }
                    idata_(i,j,k,ikhi) = -1;
                }

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    if (idata_(i,j,k,ikus)==2) {
                        if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, phi, ikus, ikhi);
                        if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, phi, ikus, ikhi);
                        if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, phi, ikus, ikhi);
                        if (i < maxi-1)
                            UpdateNeighbor(theHeap, i+1, j, k, phi, ikus, ikhi);
                        if (j < maxj-1)
                            UpdateNeighbor(theHeap, i, j+1, k, phi, ikus, ikhi);
                        if (k < maxk-1)
                            UpdateNeighbor(theHeap, i, j, k+1, phi, ikus, ikhi);
                    }
                }

        {
            int i, j, k;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,k,ikhi);
                idata_(i,j,k,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, phi, ikus, ikhi);
                if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, phi, ikus, ikhi);
                if (i < maxi-1)
                    UpdateNeighbor(theHeap, i+1, j, k, phi, ikus, ikhi);
                if (j < maxj-1)
                    UpdateNeighbor(theHeap, i, j+1, k, phi, ikus, ikhi);
                if (k < maxk-1)
                    UpdateNeighbor(theHeap, i, j, k+1, phi, ikus, ikhi);
            }
        }

        switch (dir) {
        case 0:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    for (k=0; k<maxk; ++k)
                        data_(i,j,k,func) = data_(i,j,k,phi);
            break;
        case 1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    for (k=0; k<maxk; ++k)
                        if (data_(i,j,k,phi) > 0) data_(i,j,k,func) = data_(i,j,k,phi);
            break;
        case -1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    for (k=0; k<maxk; ++k)
                        if (data_(i,j,k,phi) < 0) data_(i,j,k,func) = data_(i,j,k,phi);
            break;
        }

        bc->Apply(func);

    }
#endif

    void UniformMesh3D::ExtendDistance(const int phi, const int ikus,
                                       const int ikhi, const int imask)
    {
        NumericTrait<double> t;

        // Initialize the accepted and tentative points

        // Initialize the ikus

        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                for (k=0; k<maxk; ++k)
                    if (idata_(i,j,k,imask)) {
                        idata_(i,j,k,ikus) = 2;
                    } else {
                        idata_(i,j,k,ikus) = 0;
                        data_(i,j,k,phi) = t.max_value;
                    }   

        // Now do fast marching to complete extension

        UniformMesh3D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                for (k=0; k<maxk; ++k) {
                    if (idata_(i,j,k,ikus)==2) {
                        if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, phi, ikus, ikhi);
                        if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, phi, ikus, ikhi);
                        if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, phi, ikus, ikhi);
                        if (i < maxi-1) 
                            UpdateNeighbor(theHeap, i+1, j, k, phi, ikus, ikhi);
                        if (j < maxj-1) 
                            UpdateNeighbor(theHeap, i, j+1, k, phi, ikus, ikhi);
                        if (k < maxk-1) 
                            UpdateNeighbor(theHeap, i, j, k+1, phi, ikus, ikhi);
                    }
                }

        {
            int i, j, k;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,k,ikhi);
                idata_(i,j,k,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, k, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, k, phi, ikus, ikhi);
                if (k > 0) UpdateNeighbor(theHeap, i, j, k-1, phi, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, k, phi, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, k, phi, ikus, ikhi);
                if (k < maxk-1) UpdateNeighbor(theHeap, i, j, k+1, phi, ikus, ikhi);
            }
        }
    }

    double UniformMesh3D::TValFormula(const int sx0, const int sx1,
                                      const int sy0, const int sy1,
                                      const int sz0, const int sz1,
                                      const double xy0, 
                                      const double x1, const double x2,
                                      const double y1, const double y2,
                                      const double z1, const double z2)
    {
        double a, b, c, d;
        double ans;
        int sgn = sign(xy0);
   
        a = 0.;
        b = 0.;
        c = 0.;
        if (sx0) {
            if (sx1) {
                a += sqr(1.5/dx);
                b += (-6*x1+1.5*x2)/dx/dx;
                c += sqr((2*x1-x2/2)/dx);
            }
            else {
                a += sqr(1/dx);
                b += -2*x1/dx/dx;
                c += sqr(x1/dx);
            }
        }
        if (sy0) {
            if (sy1) {
                a += sqr(1.5/dy);
                b += (-6*y1+1.5*y2)/dy/dy;
                c += sqr((2*y1-y2/2)/dy);
            }
            else {
                a += sqr(1/dy);
                b += -2*y1/dy/dy;
                c += sqr(y1/dy);
            }
        }
        if (sz0) {
            if (sz1) {
                a += sqr(1.5/dz);
                b += (-6*z1+1.5*z2)/dz/dz;
                c += sqr((2*z1-z2/2)/dz);
            }
            else {
                a += sqr(1/dz);
                b += -2*z1/dz/dz;
                c += sqr(z1/dz);
            }
        }
        c -= 1;
        d = b*b-4*a*c;
        if (d >= 0.) {
            ans = (-b+sgn*sqrt(d))/2/a;
        }
        else {
            if (sx0 == 0) {
                if (fabs(y1) < fabs(z1)) {
                    ans = TValFormula(0,0,sy0,sy1,0,0,xy0,x1,x2,y1,y2,z1,z2);
                }
                else {
                    ans = TValFormula(0,0,0,0,sz0,sz1,xy0,x1,x2,y1,y2,z1,z2);
                }
            }
            else if (sy0 == 0) {
                if (fabs(x1) < fabs(z1)) {
                    ans = TValFormula(sx0,sx1,0,0,0,0,xy0,x1,x2,y1,y2,z1,z2);
                }
                else {
                    ans = TValFormula(0,0,0,0,sz0,sz1,xy0,x1,x2,y1,y2,z1,z2);
                }
            }
            else if (sz0 == 0) {
                if (fabs(x1) < fabs(y1)) {
                    ans = TValFormula(sx0,sx1,0,0,0,0,xy0,x1,x2,y1,y2,z1,z2);
                }
                else {
                    ans = TValFormula(0,0,sy0,sy1,0,0,xy0,x1,x2,y1,y2,z1,z2);
                }
            }
            else {
                if (fabs(x1) < fabs(z1)) {
                    if (fabs(y1) < fabs(z1)) {
                        ans = TValFormula(sx0,sx1,sy0,sy1,0,0,xy0,x1,x2,y1,y2,z1,z2);
                    }
                    else {
                        ans = TValFormula(sx0,sx1,0,0,sz0,sz1,xy0,x1,x2,y1,y2,z1,z2);
                    }
                }
                else {
                    if (fabs(y1) < fabs(x1)) {
                        ans = TValFormula(0,0,sy0,sy1,sz0,sz1,xy0,x1,x2,y1,y2,z1,z2);
                    }
                    else {
                        ans = TValFormula(sx0,sx1,0,0,sz0,sz1,xy0,x1,x2,y1,y2,z1,z2);
                    }
                }
            }
        }
   
        return ans;
    }

    double UniformMesh3D::ComputeTVal(const int i, const int j, const int k,
                                      const int l, const int ikus)
    {
        int sx[2][2];
        int sy[2][2];
        int sz[2][2];

        sx[0][0] = idata_(i-1,j,k,ikus) == 2 ? 1 : 0;
        sx[0][1] = (idata_(i-2,j,k,ikus) == 2 && sx[0][0]) ? 1 : 0;
        sx[1][0] = idata_(i+1,j,k,ikus) == 2 ? 1 : 0;
        sx[1][1] = (idata_(i+2,j,k,ikus) == 2 && sx[1][0]) ? 1 : 0;
        sy[0][0] = idata_(i,j-1,k,ikus) == 2 ? 1 : 0;
        sy[0][1] = (idata_(i,j-2,k,ikus) == 2 && sy[0][0]) ? 1 : 0;
        sy[1][0] = idata_(i,j+1,k,ikus) == 2 ? 1 : 0;
        sy[1][1] = (idata_(i,j+2,k,ikus) == 2 && sy[1][0]) ? 1 : 0;
        sz[0][0] = idata_(i,j,k-1,ikus) == 2 ? 1 : 0;
        sz[0][1] = (idata_(i,j,k-2,ikus) == 2 && sz[0][0]) ? 1 : 0;
        sz[1][0] = idata_(i,j,k+1,ikus) == 2 ? 1 : 0;
        sz[1][1] = (idata_(i,j,k+2,ikus) == 2 && sz[1][0]) ? 1 : 0;
   
        double temptval;
        double tval = DBL_MAX;

        if (sx[0][0] || sy[0][0] || sz[0][0]) {
            temptval = TValFormula(sx[0][0],sx[0][1],sy[0][0],sy[0][1],
                                   sz[0][0],sz[0][1],data_(i,j,k,l),data_(i-1,j,k,l),
                                   data_(i-2,j,k,l),data_(i,j-1,k,l),
                                   data_(i,j-2,k,l),data_(i,j,k-1,l),
                                   data_(i,j,k-2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[1][0] || sy[0][0] || sz[0][0]) {
            temptval = TValFormula(sx[1][0],sx[1][1],sy[0][0],sy[0][1],
                                   sz[0][0],sz[0][1],data_(i,j,k,l),data_(i+1,j,k,l),
                                   data_(i+2,j,k,l),data_(i,j-1,k,l),
                                   data_(i,j-2,k,l),data_(i,j,k-1,l),
                                   data_(i,j,k-2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[0][0] || sy[1][0] || sz[0][0]) {
            temptval = TValFormula(sx[0][0],sx[0][1],sy[1][0],sy[1][1],
                                   sz[0][0],sz[0][1],data_(i,j,k,l),data_(i-1,j,k,l),
                                   data_(i-2,j,k,l),data_(i,j+1,k,l),
                                   data_(i,j+2,k,l),data_(i,j,k-1,l),
                                   data_(i,j,k-2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[1][0] || sy[1][0] || sz[0][0]) {
            temptval = TValFormula(sx[1][0],sx[1][1],sy[1][0],sy[1][1],
                                   sz[0][0],sz[0][1],data_(i,j,k,l),data_(i+1,j,k,l),
                                   data_(i+2,j,k,l),data_(i,j+1,k,l),
                                   data_(i,j+2,k,l),data_(i,j,k-1,l),
                                   data_(i,j,k-2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[0][0] || sy[0][0] || sz[1][0]) {
            temptval = TValFormula(sx[0][0],sx[0][1],sy[0][0],sy[0][1],
                                   sz[1][0],sz[1][1],data_(i,j,k,l),data_(i-1,j,k,l),
                                   data_(i-2,j,k,l),data_(i,j-1,k,l),
                                   data_(i,j-2,k,l),data_(i,j,k+1,l),
                                   data_(i,j,k+2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[1][0] || sy[0][0] || sz[1][0]) {
            temptval = TValFormula(sx[1][0],sx[1][1],sy[0][0],sy[0][1],
                                   sz[1][0],sz[1][1],data_(i,j,k,l),data_(i+1,j,k,l),
                                   data_(i+2,j,k,l),data_(i,j-1,k,l),
                                   data_(i,j-2,k,l),data_(i,j,k+1,l),
                                   data_(i,j,k+2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[0][0] || sy[1][0] || sz[1][0]) {
            temptval = TValFormula(sx[0][0],sx[0][1],sy[1][0],sy[1][1],
                                   sz[1][0],sz[1][1],data_(i,j,k,l),data_(i-1,j,k,l),
                                   data_(i-2,j,k,l),data_(i,j+1,k,l),
                                   data_(i,j+2,k,l),data_(i,j,k+1,l),
                                   data_(i,j,k+2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
        if (sx[1][0] || sy[1][0] || sz[1][0]) {
            temptval = TValFormula(sx[1][0],sx[1][1],sy[1][0],sy[1][1],
                                   sz[1][0],sz[1][1],data_(i,j,k,l),data_(i+1,j,k,l),
                                   data_(i+2,j,k,l),data_(i,j+1,k,l),
                                   data_(i,j+2,k,l),data_(i,j,k+1,l),
                                   data_(i,j,k+2,l));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
      
        return data_(i,j,k,l) = tval;
    }

   
}




