/************************************************************
    $Header: um2advect.cpp,v 2.3 99/01/06 13:58:58 chopp Exp $
 
    $Log:       um2advect.cpp,v $
    Revision 2.3  99/01/06  13:58:58  13:58:58  chopp (David Chopp)
    *** none ***

    Revision 1.7  98/05/18  14:06:12  14:06:12  chopp (David Chopp)
    *** empty log message ***
  
    **********************************************************/
#include "numtrait.h"
#include "datavec.h"
#include "uniformmesh2d.h"
#include "um2heap.h"
#include "debug.h"
#include "plotwindow2d.h"
#include "um2boundary.h"
#include "utility.h"

namespace levelset {
        
    double DistToSeg(const Segment& s, const double x, const double y);

    void UniformMesh2D::Advance(const int k, const int kv, const int knorm,
                                const double dt)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,k) += dt*data_(i,j,knorm)*data_(i,j,kv);
        bc->Apply(k);
    }

    void UniformMesh2D::Advance(const int k, const int kv, const double dt)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,k) += dt*data_(i,j,kv);
        bc->Apply(k);
    }

    void UniformMesh2D::Product(const int k, const int ka, const int kb)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,k) = data_(i,j,ka)*data_(i,j,kb);
        bc->Apply(k);
    }


    double UniformMesh2D::ComputeTVal(const int i, const int j, const int k)
    {
        double a, b;
        if (data_(i,j,k) > 0) {
            if (i <= 0)
                a = data_(i+1,j,k);
            else if (i >= maxi-1)
                a = data_(i-1,j,k);
            else
                a = min(data_(i-1,j,k),data_(i+1,j,k));

            if (j <= 0)
                b = data_(i,j+1,k);
            else if (j >= maxj-1)
                b = data_(i,j-1,k);
            else
                b = min(data_(i,j-1,k),data_(i,j+1,k));

            if (b > a) {
                if (b-a < dx) {
                    // DOUBLE CHECK THESE FORMULAE
                    data_(i,j,k) = (a*dx*dx+b*dy*dy+dx*dy*sqrt(dx*dx+dy*dy-sqr(b-a)))
                        /(dx*dx+dy*dy);
                }
                else {
                    data_(i,j,k) = a+dx;
                }
            }
            else {
                if (a-b < dy) {
                    data_(i,j,k) = (a*dx*dx+b*dy*dy+dx*dy*sqrt(dx*dx+dy*dy-sqr(b-a)))
                        /(dx*dx+dy*dy);
                }
                else {
                    data_(i,j,k) = b+dy;
                }
            }
        }
        else {
            if (i > 0) {
                if (i < maxi-1) {
                    a = max(data_(i-1,j,k),data_(i+1,j,k));
                }
                else {
                    a = data_(i-1,j,k);
                }
            }
            else {
                a = data_(i+1,j,k);
            }
   
            if (j > 0) {
                if (j < maxj-1) {
                    b = max(data_(i,j-1,k),data_(i,j+1,k));
                }
                else {
                    b = data_(i,j-1,k);
                }
            }
            else {
                b = data_(i,j+1,k);
            }

            if (b < a) {
                if (a-b < dx) {
                    data_(i,j,k) = (a*dx*dx+b*dy*dy-dx*dy*sqrt(dx*dx+dy*dy-sqr(b-a)))
                        /(dx*dx+dy*dy);
                }
                else {
                    data_(i,j,k) = a-dx;
                }
            }
            else {
                if (b-a < dy) {
                    data_(i,j,k) = (a*dx*dx+b*dy*dy-dx*dy*sqrt(dx*dx+dy*dy-sqr(b-a)))
                        /(dx*dx+dy*dy);
                }
                else {
                    data_(i,j,k) = b-dy;
                }
            }
        }
        return data_(i,j,k);
    }


    void UniformMesh2D::ExtendVals(const int i, const int j, const int k,
                                   const int ikus, const int ikhi,
                                   const int* p, const int pl)
    {
        int d;
        //int sgn = data_(i,j,k) >= 0. ? 1 : -1;
        d  = ((i > 0) && idata_(i-1,j,ikus)==2) ? 1 : 0;
        d += ((j > 0) && idata_(i,j-1,ikus)==2) ? 2 : 0;
        d += ((i < maxi-1) && idata_(i+1,j,ikus)==2) ? 4 : 0;
        d += ((j < maxj-1) && idata_(i,j+1,ikus)==2) ? 8 : 0;
   
        double a;
        double b;
        int l;
        switch(d) {
        case 0:
            std::cout << "Case 0 hit\n";
            a = (data_(i+1,j,k)+data_(i,j+1,k)+data_(i-1,j,k)+data_(i,j-1,k)
                 -4.*data_(i,j,k))/4.;
            a = 0.;
            b = 0.;
            for (l=0; l<pl; ++l) {
                if (idata_(i+1,j,ikhi)) {
                    a += data_(i+1,j,p[l]);
                    b += 1.;
                }
                if (idata_(i-1,j,ikhi)) {
                    a += data_(i-1,j,p[l]);
                    b += 1.;
                }
                if (idata_(i,j+1,ikhi)) {
                    a += data_(i,j+1,p[l]);
                    b += 1.;
                }
                if (idata_(i,j-1,ikhi)) {
                    a += data_(i,j-1,p[l]);
                    b += 1.;
                }
                data_(i,j,p[l]) = a/b;
            }
            break;
        case 1:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i-1,j,p[l]);
            break;
        case 2:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i,j-1,p[l]);
            break;
        case 3:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)
                                          -data_(i,j-1,k))*data_(i,j-1,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(i-1,j,k))
                                   *data_(i-1,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                      +dy*dy*(data_(i,j,k)-data_(i-1,j,k)));
            break;
        case 4:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i+1,j,p[l]);
            break;
        case 5:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i-1,j,p[l])*(data_(i,j,k)-data_(i-1,j,k))
                                   -data_(i+1,j,p[l])*(data_(i+1,j,k)-data_(i,j,k)))
                    /(2.*data_(i,j,k)-data_(i-1,j,k)-data_(i+1,j,k));
            break;
        case 6:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                                   *data_(i,j-1,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(i+1,j,k))
                                   *data_(i+1,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                      +dy*dy*(data_(i,j,k)-data_(i+1,j,k)));
            break;
        case 7:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                                   *data_(i,j-1,p[l])
                                   +dy*dy*(data_(i-1,j,p[l])
                                           *(data_(i,j,k)-data_(i-1,j,k))
                                           -data_(i+1,j,p[l])
                                           *(data_(i+1,j,k)-data_(i,j,k))))
                    /(dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                      +dy*dy*(2*data_(i,j,k)-data_(i-1,j,k)
                              -data_(i+1,j,k)));
            break;
        case 8:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i,j+1,p[l]);
            break;
        case 9:
            for (l=0; l<pl; ++l) {
                a = dx*dx*(data_(i,j,k)-data_(i,j+1,k));
                b = dy*dy*(data_(i,j,k)-data_(i-1,j,k));
                data_(i,j,p[l]) = (a*data_(i,j+1,p[l])+b*data_(i-1,j,p[l]))/(a+b);
            }
            break;
        case 10:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i,j-1,p[l])*(data_(i,j,k)-data_(i,j-1,k))
                                   -data_(i,j+1,p[l])*(data_(i,j+1,k)
                                                       -data_(i,j,k)))
                    /(2*data_(i,j,k)-data_(i,j-1,k)-data_(i,j+1,k));
            break;
        case 11:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j+1,p[l])
                                          *(data_(i,j,k)-data_(i,j+1,k))
                                          +data_(i,j-1,p[l])
                                          *(data_(i,j,k)-data_(i,j-1,k)))
                                   +dy*dy*data_(i-1,j,p[l])
                                   *(data_(i,j,k)-data_(i-1,j,k)))
                    /(dx*dx*(2*data_(i,j,k)
                             -data_(i,j-1,k)-data_(i,j+1,k))
                      +dy*dy*(data_(i,j,k)-data_(i-1,j,k)));
            break;
        case 12:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)
                                          -data_(i,j+1,k))*data_(i,j+1,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(i+1,j,k))
                                   *data_(i+1,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,j+1,k))
                      +dy*dy*(data_(i,j,k)-data_(i+1,j,k)));
            break;
        case 13:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*data_(i,j+1,p[l])
                                   *(data_(i,j,k)-data_(i,j+1,k))
                                   +dy*dy*(data_(i-1,j,p[l])
                                           *(data_(i,j,k)-data_(i-1,j,k))
                                           +data_(i+1,j,p[l])
                                           *(data_(i,j,k)-data_(i+1,j,k))))
                    /(dx*dx*(data_(i,j,k)-data_(i,j+1,k))
                      +dy*dy*(2*data_(i,j,k)
                              -data_(i-1,j,k)-data_(i+1,j,k)));
            break;
        case 14:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j-1,p[l])
                                          *(data_(i,j,k)-data_(i,j-1,k))
                                          +data_(i,j+1,p[l])
                                          *(data_(i,j,k)-data_(i,j+1,k)))
                                   +dy*dy*data_(i+1,j,p[l])
                                   *(data_(i,j,k)-data_(i+1,j,k)))
                    /(dx*dx*(2*data_(i,j,k)-data_(i,j-1,k)
                             -data_(i,j+1,k))
                      +dy*dy*(data_(i,j,k)-data_(i+1,j,k)));
            break;
        case 15:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i-1,j,p[l])+data_(i,j-1,p[l])
                                   +data_(i+1,j,p[l])+data_(i,j+1,p[l]))/4.;
            break;
        }
    }


    void UniformMesh2D::ExtendFunc(const int i, const int j, const int k,
                                   const int ikus)
    {
        int d;
        double temp[4];
        int sgn = data_(i,j,k) >= 0. ? 1 : -1;
        d  = ((i > 0) && (idata_(i-1,j,ikus) == 2)) ? 1 : 0;
        d += ((j > 0) && (idata_(i,j-1,ikus) == 2)) ? 2 : 0;
        d += ((i < maxi-1) && (idata_(i+1,j,ikus) == 2)) ? 4 : 0;
        d += ((j < maxj-1) && (idata_(i,j+1,ikus) == 2)) ? 8 : 0;
   
        switch(d) {
        case 0:
            break;
        case 1:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                data_(i,j,k) = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                data_(i,j,k) = sgn*dx+data_(i-1,j,k);
            break;
        case 2:
            if (j>1 && idata_(i,j-2,ikus) == 2)
                data_(i,j,k) = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                data_(i,j,k) = sgn*dy+data_(i,j-1,k);
            break;
        case 3:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[1] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[1] = sgn*dy+data_(i,j-1,k);
            data_(i,j,k) = (dy*dy*temp[0]+dx*dx*temp[1])/(dx*dx+dy*dy);
            break;
        case 4:
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                data_(i,j,k) = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                data_(i,j,k) = sgn*dx+data_(i+1,j,k);
            break;
        case 5:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[1] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[1] = sgn*dx+data_(i+1,j,k);
            data_(i,j,k) = (temp[0]+temp[1])/2.;
            break;
        case 6:
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[0] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[0] = sgn*dy+data_(i,j-1,k);
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[1] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[1] = sgn*dx+data_(i+1,j,k);
            data_(i,j,k) = (dy*dy*temp[1]+dx*dx*temp[0])/(dx*dx+dy*dy);
            break;
        case 7:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[1] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[1] = sgn*dy+data_(i,j-1,k);
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[2] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[2] = sgn*dx+data_(i+1,j,k);
            data_(i,j,k) = (dy*dy*temp[0]+dx*dx*temp[1]+dy*dy*temp[2])
                /(dx*dx+2.*dy*dy);
            break;
        case 8:
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                data_(i,j,k) = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                data_(i,j,k) = sgn*dy+data_(i,j+1,k);
            break;
        case 9:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[1] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[1] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (dy*dy*temp[0]+dx*dx*temp[1])/(dx*dx+dy*dy);
            break;
        case 10:
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[0] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[0] = sgn*dy+data_(i,j-1,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[1] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[1] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (temp[0]+temp[1])/2.;
            break;
        case 11:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[1] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[1] = sgn*dy+data_(i,j-1,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[2] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[2] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (dy*dy*temp[0]+dx*dx*temp[1]+dx*dx*temp[2])
                /(dx*dx+2.*dy*dy);
            break;
        case 12:
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[0] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[0] = sgn*dx+data_(i+1,j,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[1] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[1] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (dy*dy*temp[0]+dx*dx*temp[1])/(dx*dx+dy*dy);
            break;
        case 13:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[1] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[1] = sgn*dx+data_(i+1,j,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[2] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[2] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (dy*dy*temp[0]+dy*dy*temp[1]+dx*dx*temp[2])
                /(dx*dx+2.*dy*dy);
            break;
        case 14:
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[0] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[0] = sgn*dy+data_(i,j-1,k);
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[1] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[1] = sgn*dx+data_(i+1,j,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[2] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[2] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (dx*dx*temp[0]+dy*dy*temp[1]+dx*dx*temp[2])
                /(dx*dx+2.*dy*dy);
            break;
        case 15:
            if (i>1 && idata_(i-2,j,ikus) == 2)
                temp[0] = 2*data_(i-1,j,k)-data_(i-2,j,k);
            else
                temp[0] = sgn*dx+data_(i-1,j,k);
            if (j>1 && idata_(i,j-2,ikus) == 2)
                temp[1] = 2*data_(i,j-1,k)-data_(i,j-2,k);
            else
                temp[1] = sgn*dy+data_(i,j-1,k);
            if (i<maxi-2 && idata_(i+2,j,ikus) == 2)
                temp[2] = 2*data_(i+1,j,k)-data_(i+2,j,k);
            else
                temp[2] = sgn*dx+data_(i+1,j,k);
            if (j<maxj-2 && idata_(i,j+2,ikus) == 2)
                temp[3] = 2*data_(i,j+1,k)-data_(i,j+2,k);
            else
                temp[3] = sgn*dy+data_(i,j+1,k);
            data_(i,j,k) = (dy*dy*temp[0]+dx*dx*temp[1]+dy*dy*temp[2]
                            +dx*dx*temp[3])
                /(dx*dx+dy*dy)/2.;
            break;
        }
    }

#if 0
    void UniformMesh2D::UpdateNeighbor(UniformMesh2D_Heap& heap, const int i,
                                       const int j, const int k, const int ikus,
                                       const int ikhi)
    {
        switch(idata_(i,j,ikus)) {
        case 0:
            heap.Insert(i,j,abs(ComputeTVal(i,j,k)),ikhi);
            idata_(i,j,ikus) = 1;
            break;
        case 1:
            heap.Change(idata_(i,j,ikhi),abs(ComputeTVal(i,j,k)),ikhi);
            break;
        case 2:
            break;
        default:
            cerr << "Error: Gridpoint " << i << "," << j << " not initialized.\n";
            break;
        }
    }
#else
    void UniformMesh2D::UpdateNeighbor(UniformMesh2D_Heap& heap, const int i,
                                       const int j, const int k, const int ikus,
                                       const int ikhi)
    {
#ifdef LEVEL_DEBUG
        int ii;
        double d;
#endif
        switch(idata_(i,j,ikus)) {
        case 0:
            heap.Insert(i,j,fabs(ComputeTVal(i,j,k,ikus)),ikhi);
            idata_(i,j,ikus) = 1;
            break;
        case 1:
#ifdef LEVEL_DEBUG
            ii = idata_(i,j,ikhi);
            d = fabs(ComputeTVal(i,j,k,ikus));
            heap.Change(ii,d,ikhi);
#else
            heap.Change(idata_(i,j,ikhi),fabs(ComputeTVal(i,j,k,ikus)),ikhi);
#endif
            break;
        case 2:
            break;
        default:
            std::cerr << "Error: Gridpoint " << i << "," << j << " not initialized.\n";
            break;
        }
    }
#endif

    void UniformMesh2D::AddInitVelocity(const int i, const int j, const int ikus,
                                        const int ikhi, const int func,
                                        const int vel, const Direction d,
                                        const double f)
    {
        if (idata_(i,j,ikhi) & d) {

            double s1, s2, t1, t2;

            if (idata_(i,j,ikus) & DLeft) 
                s1 = frac(0., data_(i,j,func), data_(i-1,j,func));
            if (idata_(i,j,ikus) & DDown) 
                t1 = frac(0., data_(i,j,func), data_(i,j-1,func));
            if (idata_(i,j,ikus) & DRight) 
                s2 = frac(0., data_(i,j,func), data_(i+1,j,func));
            if (idata_(i,j,ikus) & DUp)
                t2 = frac(0., data_(i,j,func), data_(i,j+1,func));
            int sgn = (data_(i,j,func) >= 0) ? 1 : -1;
        
                        
            switch(d) {
            case DLeft:
                switch(idata_(i,j,ikus)) {
                case DLeft: 
                    data_(i,j,vel) += f;
                    break;
                case DLeft | DDown:
                    data_(i,j,vel) += sqr(t1*dy)*f/(sqr(s1*dx)+sqr(t1*dy));
                    break;
                case DLeft | DRight:
                    if (!(idata_(i,j,ikhi) & DRight)) {
                        double f2 = data_(i,j,vel);
                        if (s2 >= s1) {
                            if (dx*(s2-s1)+sgn*(f2-f) >= 0) {
                                data_(i,j,vel) = f;
                            } else {
                                data_(i,j,vel) = f+dx*(s2-s1);
                            }
                        } else {
                            if (dx*(s2-s1)+sgn*(f2-f) >= 0) {
                                data_(i,j,vel) = f2-dx*(s2-s1);
                            } else {
                                data_(i,j,vel) = f2;
                            }
                        }
                    } else {
                        data_(i,j,vel) = f;
                    }
                    break;
                case DLeft | DDown | DRight:
                    data_(i,j,vel) += sqr(t1*dy)*s2*(s2-s1)*f
                        /(sqr(s1*s2*dx)+sqr(t1*(s2-s1)*dy));
                    break;
                case DLeft | DUp:
                    data_(i,j,vel) += sqr(t2*dy)*f/(sqr(s1*dx)+sqr(t2*dy));
                    break;
                case DLeft | DDown | DUp:
                    data_(i,j,vel) += sqr(t1*t2*dy)*f
                        /(sqr(s1*(t2-t1)*dx)+sqr(t1*t2*dy));
                    break;
                case DLeft | DRight | DUp:
                    data_(i,j,vel) += sqr(t2*dy)*s2*(s2-s1)*f
                        /(sqr(s1*s2*dx)+sqr(t2*(s2-s1)*dy));
                    break;
                case DUp | DLeft | DDown | DRight:
                    data_(i,j,vel) += sqr(t1*t2*dx)*s2*(s2-s1)*f
                        /(sqr(s1*s2*(t2-t1)*dx)+sqr(t1*t2*(s2-s1)*dy));
                    break;
                }
                break;
            case DDown:
                switch(idata_(i,j,ikus)) {
                case DDown: 
                    data_(i,j,vel) += f;
                    break;
                case DDown | DLeft: 
                    data_(i,j,vel) += sqr(s1*dx)*f/(sqr(s1*dx)+sqr(t1*dy));
                    break;
                case DDown | DRight:
                    data_(i,j,vel) += sqr(s2*dx)*f/(sqr(s2*dx)+sqr(t1*dy));
                    break;
                case DDown | DLeft | DRight:
                    data_(i,j,vel) += sqr(s1*s2*dx)*f
                        /(sqr(s1*s2*dx)+sqr(t1*(s2-s1)*dy));
                    break;
                case DDown | DUp:
                    if (!(idata_(i,j,ikhi) & DUp)) {
                        double f2 = data_(i,j,vel);
                        if (t2 >= t1) {
                            if (dy*(t2-t1)+sgn*(f2-f) >= 0) {
                                data_(i,j,vel) = f;
                            } else {
                                data_(i,j,vel) = f+dy*(t2-t1);
                            }
                        } else {
                            if (dy*(t2-t1)+sgn*(f2-f) >= 0) {
                                data_(i,j,vel) = f2-dy*(t2-t1);
                            } else {
                                data_(i,j,vel) = f2;
                            }
                        }
                    } else {
                        data_(i,j,vel) = f;
                    }
                    break;
                case DDown | DUp | DLeft:
                    data_(i,j,vel) += sqr(s1*dx)*t2*(t2-t1)*f
                        /(sqr(s1*(t2-t1)*dx)+sqr(t1*t2*dy));
                    break;
                case DDown | DUp | DRight:
                    data_(i,j,vel) += sqr(s2*dx)*t2*(t2-t1)*f
                        /(sqr(s2*(t2-t1)*dx)+sqr(t1*t2*dy));
                    break;
                case DUp | DLeft | DDown | DRight:
                    data_(i,j,vel) += sqr(s1*s2*dx)*t2*(t2-t1)*f
                        /(sqr(s1*s2*(t2-t1)*dx)+sqr(t1*t2*(s2-s1)*dy));
                    break;
                }
                break;
            case DRight:
                switch(idata_(i,j,ikus)) {
                case DRight: 
                    data_(i,j,vel) += f;
                    break;
                case DRight | DLeft:
                    if (!(idata_(i,j,ikhi) & DLeft)) {
                        double f2 = data_(i,j,vel);
                        if (s2 >= s1) {
                            if (dx*(s2-s1)+sgn*(f-f2) >= 0) {
                                data_(i,j,vel) = f2;
                            } else {
                                data_(i,j,vel) = f2+dx*(s2-s1);
                            }
                        } else {
                            if (dx*(s2-s1)+sgn*(f-f2) >= 0) {
                                data_(i,j,vel) = f-dx*(s2-s1);
                            } else {
                                data_(i,j,vel) = f;
                            }
                        }
                    } else {
                        data_(i,j,vel) = f;
                    }
                    break;
                case DRight | DDown:
                    data_(i,j,vel) += sqr(t1*dy)*f/(sqr(s2*dx)+sqr(t1*dy));
                    break;
                case DRight | DLeft | DDown:
                    data_(i,j,vel) -= sqr(t1*dy)*s1*(s2-s1)*f
                        /(sqr(s1*s2*dx)+sqr(t1*(s2-s1)*dy));
                    break;
                case DRight | DUp:
                    data_(i,j,vel) += sqr(t2*dy)*f/(sqr(s2*dx)+sqr(t2*dy));
                    break;
                case DRight | DLeft | DUp:
                    data_(i,j,vel) -= sqr(t2*dy)*s1*(s2-s1)*f
                        /(sqr(s1*s2*dx)+sqr(t2*(s2-s1)*dy));
                    break;
                case DRight | DDown | DUp:
                    data_(i,j,vel) += sqr(t1*t2*dy)*f
                        /(sqr(s2*(t2-t1)*dx)+sqr(t1*t2*dy));
                    break;
                case DUp | DLeft | DDown | DRight:
                    data_(i,j,vel) -= sqr(t1*t2*dy)*s1*(s2-s1)*f
                        /(sqr(s1*s2*(t2-t1)*dx)+sqr(t1*t2*(s2-s1)*dy));
                    break;
                }
                break;
            case DUp:
                switch(idata_(i,j,ikus)) {
                case DUp: 
                    data_(i,j,vel) += f;
                    break;
                case DUp | DLeft:
                    data_(i,j,vel) += sqr(s1*dx)*f/(sqr(s1*dx)+sqr(t2*dy));
                    break;
                case DUp | DDown:
                    if (!(idata_(i,j,ikhi) & DDown)) {
                        double f2 = data_(i,j,vel);
                        if (t2 >= t1) {
                            if (dy*(t2-t1)+sgn*(f-f2) >= 0) {
                                data_(i,j,vel) = f2;
                            } else {
                                data_(i,j,vel) = f2+dy*(t2-t1);
                            }
                        } else {
                            if (dy*(t2-t1)+sgn*(f-f2) >= 0) {
                                data_(i,j,vel) = f-dy*(t2-t1);
                            } else {
                                data_(i,j,vel) = f;
                            }
                        }
                    } else {
                        data_(i,j,vel) = f;
                    }
                    break;
                case DUp | DLeft | DDown:
                    data_(i,j,vel) -= sqr(s1*dx)*t1*(t2-t1)*f
                        /(sqr(s1*(t2-t1)*dx)+sqr(t1*t2*dy));
                    break;
                case DUp | DRight:
                    data_(i,j,vel) += sqr(s2*dx)*f/(sqr(s2*dx)+sqr(t2*dy));
                    break;
                case DUp | DLeft | DRight:
                    data_(i,j,vel) += sqr(s1*s2*dx)*f
                        /(sqr(s1*s2*dx)+sqr(t2*(s2-s1)*dy));
                    break;
                case DUp | DDown | DRight:
                    data_(i,j,vel) -= sqr(s2*dx)*t1*(t2-t1)*f
                        /(sqr(s2*(t2-t1)*dx)+sqr(t1*t2*dy));
                    break;
                case DUp | DLeft | DDown | DRight:
                    data_(i,j,vel) -= sqr(s1*s2*dx)*t1*(t2-t1)*f
                        /(sqr(s1*s2*(t2-t1)*dx)+sqr(t1*t2*(s2-s1)*dy));
                    break;
                }
                break;
            }
            idata_(i,j,ikhi) -= d;
        }
    }

    void UniformMesh2D::InitLevelSet(const int i, const int j, const int func, 
                                     const int phi, const int ikus)
    {
        Bicubic p(data_(i,j,func), data_(i+1,j,func), data_(i,j+1,func), 
                  data_(i+1,j+1,func), (data_(i+1,j,func)-data_(i-1,j,func))/2/dx, 
                  (data_(i+2,j,func)-data_(i,j,func))/2/dx, 
                  (data_(i+1,j+1,func)-data_(i-1,j+1,func))/2/dx, 
                  (data_(i+2,j+1,func)-data_(i,j+1,func))/2/dx,
                  (data_(i,j+1,func)-data_(i,j-1,func))/2/dy, 
                  (data_(i+1,j+1,func)-data_(i+1,j-1,func))/2/dy,
                  (data_(i,j+2,func)-data_(i,j,func))/2/dy, 
                  (data_(i+1,j+2,func)-data_(i+1,j,func))/2/dy,
                  (data_(i+1,j+1,func)+data_(i-1,j-1,func)
                   -data_(i-1,j+1,func)-data_(i+1,j-1,func))/4/dx/dy,
                  (data_(i+2,j+1,func)+data_(i,j-1,func)-data_(i,j+1,func)
                   -data_(i+2,j-1,func))/4/dx/dy,
                  (data_(i+1,j+2,func)+data_(i-1,j,func)-data_(i-1,j+2,func)
                   -data_(i+1,j,func))/4/dx/dy,
                  (data_(i+2,j+2,func)+data_(i,j,func)-data_(i,j+2,func)
                   -data_(i+2,j,func))/4/dx/dy, dx, dy);
        
        double ax, ay;
        double dist[2][2];
        char clean[2][2];
        dist[0][0] = p.LocalDist(0., 0., ax, ay, clean[0][0]);
        dist[1][0] = p.LocalDist(dx, 0., ax, ay, clean[1][0]);
        dist[0][1] = p.LocalDist(0., dy, ax, ay, clean[0][1]);
        dist[1][1] = p.LocalDist(dx, dy, ax, ay, clean[1][1]);
   
        for (int ii=0; ii<2; ++ii)
            for (int jj=0; jj<2; ++jj) {
                if (clean[ii][jj]) {
                    if (idata_(i+ii,j+jj,ikus) < 2) {
                        data_(i+ii,j+jj,phi) = copysign(dist[ii][jj],
                                                        data_(i+ii,j+jj,func));
                        idata_(i+ii,j+jj,ikus) = 2;
                    }
                    else if (fabs(data_(i+ii,j+jj,phi)) > dist[ii][jj]) {
                        data_(i+ii,j+jj,phi) =copysign(dist[ii][jj],
                                                       data_(i+ii,j+jj,func));
                        idata_(i+ii,j+jj,ikus) = 2;
                    }
                } else if (fabs(data_(i+ii,j+jj,phi)) > dist[ii][jj]) {
                    data_(i+ii,j+jj,phi) =copysign(dist[ii][jj],data_(i+ii,j+jj,func));
                    idata_(i+ii,j+jj,ikus) = 1;
                }
            }
    }   

    void UniformMesh2D::InitLevelSet(const int i, const int j, const int func,
                                     const int phi, const int ikus, const int ikhi)
    {
        if (idata_(i,j,ikhi) == idata_(i,j,ikus)) {

            double s1, s2, t1, t2;

            if (idata_(i,j,ikus) & DLeft) 
                s1 = frac(0., data_(i,j,func), data_(i-1,j,func));
            if (idata_(i,j,ikus) & DDown) 
                t1 = frac(0., data_(i,j,func), data_(i,j-1,func));
            if (idata_(i,j,ikus) & DRight) 
                s2 = frac(0., data_(i,j,func), data_(i+1,j,func));
            if (idata_(i,j,ikus) & DUp)
                t2 = frac(0., data_(i,j,func), data_(i,j+1,func));
            int sgn = (data_(i,j,func) >= 0) ? 1 : -1;
        
                        
            switch(idata_(i,j,ikus)) {
            case 0: break;
            case 1: data_(i,j,phi) = sgn*s1*dx; break;
            case 2: data_(i,j,phi) = sgn*t1*dy; break;
            case 3: data_(i,j,phi) = sgn/sqrt(1./sqr(s1*dx)+1./sqr(t1*dy)); break;
            case 4: data_(i,j,phi) = sgn*s2*dx; break;
            case 5: data_(i,j,phi) = sgn*min(s1*dx,s2*dx); break;
            case 6: data_(i,j,phi) = sgn/sqrt(1./sqr(t1*dy) +1./sqr(s2*dx)); break;
            case 7: data_(i,j,phi) = sgn*min(1./sqrt(1./sqr(t1*dy)+1./sqr(s2*dx)),
                                             1./sqrt(1./sqr(t1*dy)+1./sqr(s1*dx)));
                break;
            case 8: data_(i,j,phi) = sgn*t2*dy; break;
            case 9: data_(i,j,phi) = sgn/sqrt(1./sqr(s1*dx)+1./sqr(t2*dy)); break;
            case 10: data_(i,j,phi) = sgn*min(t1*dy,t2*dy); break;
            case 11: data_(i,j,phi) = sgn*min(1./sqrt(1./sqr(s1*dx)+1./sqr(t1*dy)),
                                              1./sqrt(1./sqr(s1*dx)+1./sqr(t2*dy)));
                break;
            case 12: data_(i,j,phi) = sgn/sqrt(1./sqr(s2*dx) +1./sqr(t2*dy)); break;
            case 13: data_(i,j,phi) = sgn*min(1./sqrt(1./sqr(s1*dx)+1./sqr(t2*dy)),
                                              1./sqrt(1./sqr(t2*dy)+1./sqr(s2*dx)));
                break;
            case 14: data_(i,j,phi) = sgn*min(1./sqrt(1./sqr(t1*dy)+1./sqr(s2*dx)),
                                              1./sqrt(1./sqr(s2*dx)+1./sqr(t2*dy)));
                break;
            case 15: data_(i,j,phi) = sgn*min(1./sqrt(1./sqr(s1*dx)+1./sqr(t1*dy)),
                                              1./sqrt(1./sqr(s1*dx)+1./sqr(t2*dy)),
                                              1./sqrt(1./sqr(s2*dx)+1./sqr(t1*dy)),
                                              1./sqrt(1./sqr(s2*dx)
                                                      +1./sqr(t2*dy)));
                break;
            }
        }
    }

    double DistToSeg(const Segment& s, const double x, const double y)
    {
        double dx[2], dy[2];
        dx[0] = s.X2()-s.X1();
        dx[1] = x-s.X1();
        dy[0] = s.Y2()-s.Y1();
        dy[1] = y-s.Y1();
        double t = (dx[0]*dx[1]+dy[0]*dy[1])/(sqr(dx[0])+sqr(dy[0]));
        if (t < 0.) 
            return sqrt(sqr(dx[1])+sqr(dy[1]));
        else if (t > 1.) 
            return sqrt(sqr(x-s.X2())+sqr(y-s.Y2()));
        else
            return fabs(-x*dy[0]+s.X1()*(s.Y2()-y)+s.X2()*dy[1])
                /(sqr(dx[0])+sqr(dy[0]));
    }

    void UniformMesh2D::ExtendVelocity(const Interface& seg, const int sl,
                                       const int sh, const int k, const int ktemp,
                                       const int ikus, const int ikhi,
                                       const int imask, const int kvel)
    {
        NumericTrait<double> t;
      
        // Initialize the accepted and tentative points
   
        // Initialize the ikus

        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                idata_(i,j,ikus) = 0;
                data_(i,j,ktemp) = data_(i,j,k) >= 0. ? t.max_value : t.min_value;
                data_(i,j,kvel) = 0.;
            }
   
        // Initialize the fast phi function and velocity values

        int ind[2];
        double mx, my;
        for (i=0; i<seg.Length(); ++i) {
            seg[i].MidPoint(mx, my);
            ind[0] = I(mx);
            ind[1] = J(my);
            if (seg[i].Left() || seg[i].Down()) 
                data_(i,j,ktemp) = sign(data_(i,j,ktemp))
                    *min(double(fabs(data_(i,j,ktemp))),
                         DistToSeg(seg[i],X(ind[0]),Y(ind[1])));
            if (seg[i].Left() || seg[i].Up()) 
                data_(i,j+1,ktemp) = sign(data_(i,j+1,ktemp))
                    *min(double(fabs(data_(i,j+1,ktemp))),
                         DistToSeg(seg[i],X(ind[0]),Y(ind[1]+1)));
            if (seg[i].Right() || seg[i].Down()) 
                data_(i+1,j,ktemp) = sign(data_(i+1,j,ktemp))
                    *min(double(fabs(data_(i+1,j,ktemp))),
                         DistToSeg(seg[i],X(ind[0]+1),Y(ind[1])));
            if (seg[i].Right() || seg[i].Up()) 
                data_(i+1,j+1,ktemp) = sign(data_(i+1,j+1,ktemp))
                    *min(double(fabs(data_(i+1,j+1,ktemp))),
                         DistToSeg(seg[i],X(ind[0]+1),Y(ind[1]+1)));
        }
   
         

        // Initialize the velocity values here
      
        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                idata_(i,j,ikhi) = -1;
            }

        for (i=-bc->Width(BDRY2D_XLO); i<maxi+bc->Width(BDRY2D_XHI); ++i)
            for (j=-bc->Width(BDRY2D_YLO); j<maxj+bc->Width(BDRY2D_YHI); ++j) {
                idata_(i,j,ikus) = 0;
                data_(i,j,ktemp) = t.max_value;
            }

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) >= 16) {
                    data_(i,j,ktemp) = 0.;
                    idata_(i,j,ikus) = 2;
                } else {
                    switch(idata_(i,j,ikus)) {
                    case 0:
                        data_(i,j,ktemp) = t.max_value;
                        break;
                    case 1:
                    case 2:
                    case 3:
                        data_(i,j,ktemp) = dx;
                        idata_(i,j,ikus) = 1;
                        theHeap.Insert(i,j,0.,ikhi);
                        break;
                    case 4:
                    case 8:
                    case 12:
                        data_(i,j,ktemp) = dy;
                        idata_(i,j,ikus) = 1;
                        theHeap.Insert(i,j,0.,ikhi);
                        break;
                    case 5:
                    case 6:
                    case 7:
                    case 9:
                    case 10:
                    case 11:
                    case 13:
                    case 14:
                    case 15:
                        data_(i,j,ktemp) = sqrt(1/dx/dx+1/dy/dy);
                        idata_(i,j,ikus) = 1;
                        theHeap.Insert(i,j,0.,ikhi);
                        break;
                    }
                }

        {
            int i, j;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask)) ExtendVals(i, j, ktemp, ikus, ikhi,
                                                   &kvel, 1);

                if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
            }
        }


    }      

    void UniformMesh2D::ExtendVelocity(const int ktemp, const int kvel,
                                       const int ikus, const int ikhi,
                                       const int imask)
    {  
        NumericTrait<double> t;
      
        // Initialize the accepted and tentative points
      
        // Initialize the UPDATE_STATUS

        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                    data_(i,j,ktemp) = 0.;
                } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,ktemp) = t.max_value;
                }
   
        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 2) {
                    if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                    if (i < maxi-1)
                        UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                    if (j < maxj-1)
                        UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                }
            
        {
            int i, j;
//         PlotWindow2D exdisp(true);
//         exdisp.SetRange(-1., maxi, -1., maxj);
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
//            exdisp.Dot(i,j);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask))
                    ExtendVals(i, j, ktemp, ikus, ikhi, &kvel, 1);

                if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
            }
        }

    }      
    
	 //overloaded for initial velocities at non-zero phi values
	 void UniformMesh2D::ExtendVelocity(const int kphi, const int ktemp, const int kvel,
                                       const int ikus, const int ikhi,
                                       const int imask)
    {  
        NumericTrait<double> t;
      
        // Initialize the accepted and tentative points
      
        // Initialize the UPDATE_STATUS

        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                    //data_(i,j,ktemp) = 0.;
                	  data_(i,j,ktemp) = fabs(data_(i,j,kphi));
					 } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,ktemp) = t.max_value;
                }
   
        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 2) {
                    if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                    if (i < maxi-1)
                        UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                    if (j < maxj-1)
                        UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                }
            
        {
            int i, j;
//         PlotWindow2D exdisp(true);
//         exdisp.SetRange(-1., maxi, -1., maxj);
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
//            exdisp.Dot(i,j);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask))
                    ExtendVals(i, j, ktemp, ikus, ikhi, &kvel, 1);

                if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
            }
        }

    }      

#if 0
    void UniformMesh2D::Reinitialize(const int func, const int phi, const int ikus,
                                     const int ikhi, const int dir)
    {
        NumericTrait<double> t;

        // Initialize the accepted and tentative points

        // Initialize the ikus

        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                idata_(i,j,ikus) = 0;
                idata_(i,j,ikhi) = 0;
                data_(i,j,phi) = data_(i,j,func) >= 0. ? t.max_value : t.min_value;
            }

        int cross;
        for (i=0; i<maxi-1; ++i)
            for (j=0; j<maxj-1; ++j) {

                if (data_(i,j,func)*data_(i+1,j,func) <= 0) {
                    idata_(i,j,ikus) += DRight;
                    idata_(i+1,j,ikus) += DLeft;
                }

                if (data_(i,j,func)*data_(i,j+1,func) <= 0) {
                    idata_(i,j,ikus) += DUp;
                    idata_(i,j+1,ikus) += DDown;
                }

                if (j==maxj-2 && data_(i,j+1,func)*data_(i+1,j+1,func) <= 0) {
                    idata_(i,j+1,ikus) += DRight;
                    idata_(i+1,j+1,ikus) += DLeft;
                }

                if (i==maxi-2 && data_(i+1,j,func)*data_(i+1,j+1,func) <= 0) {
                    idata_(i+1,j,ikus) += DUp;
                    idata_(i+1,j+1,ikus) += DDown;
                }
            }

        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                idata_(i,j,ikhi) = idata_(i,j,ikus);
                if (idata_(i,j,ikus)) {
                    InitLevelSet(i,j,func,phi,ikus,ikhi);
                    idata_(i,j,ikus) = 2;
                }
                idata_(i,j,ikhi) = -1;
            }

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (idata_(i,j,ikus)==2) {
                    if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                    if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                    if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
                }
            }

        {
            int i, j;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
                idata_(i,j,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
            }
        }

        switch (dir) {
        case 0:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    data_(i,j,func) = data_(i,j,phi);
            break;
        case 1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    if (data_(i,j,phi) > 0) data_(i,j,func) = data_(i,j,phi);
            break;
        case -1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    if (data_(i,j,phi) < 0) data_(i,j,func) = data_(i,j,phi);
            break;
        }

        bc->Apply(func);

    }

#else

    void UniformMesh2D::Reinitialize(const int func, const int phi, const int ikus,
                                     const int ikhi, const int dir)
    {
        NumericTrait<double> t;

        // Initialize the accepted and tentative points

        // Initialize the ikus

        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                idata_(i,j,ikus) = 0;
                idata_(i,j,ikhi) = -1;
                data_(i,j,phi) = data_(i,j,func) >= 0. ? t.max_value : t.min_value;
            }

        int cross;
        for (i=0; i<maxi-1; ++i)
            for (j=0; j<maxj-1; ++j) {
                cross = data_(i,j,func) > 0. ? 1 : 0;
                cross += data_(i+1,j,func) > 0 ? 2 : 0;
                cross += data_(i,j+1,func) > 0 ? 4 : 0;
                cross += data_(i+1,j+1,func) > 0 ? 8 : 0;
                if (cross > 0 && cross < 15) 
                    InitLevelSet(i,j,func,phi,ikus);
            }
        
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                if (idata_(i,j,ikus) == 1) idata_(i,j,ikus) = 0;
                
        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (idata_(i,j,ikus)==2) {
                    if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                    if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                    if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
                }
            }

        {
            int i, j;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
                idata_(i,j,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
            }
        }

        switch (dir) {
        case 0:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    data_(i,j,func) = data_(i,j,phi);
            break;
        case 1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    if (data_(i,j,phi) > 0) data_(i,j,func) = data_(i,j,phi);
            break;
        case -1:
            for (i=0; i<maxi; ++i)
                for (j=0; j<maxj; ++j)
                    if (data_(i,j,phi) < 0) data_(i,j,func) = data_(i,j,phi);
            break;
        }

        bc->Apply(func);

    }
#endif

    void UniformMesh2D::ExtendDistance(const int phi, const int ikus,
                                       const int ikhi, const int imask)
    {
        NumericTrait<double> t;

        // Initialize the accepted and tentative points

        // Initialize the ikus

        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,phi) = data_(i,j,phi) > 0. ? t.max_value : -t.max_value;
                }   

        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);

        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (idata_(i,j,ikus)==2) {
                    if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                    if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                    if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
                }
            }

        {
            int i, j;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
                idata_(i,j,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
            }
        }
    }

    double UniformMesh2D::TValFormula(const int sx0, const int sx1,
                                      const int sy0, const int sy1,
                                      const double xy0, 
                                      const double x1, const double x2,
                                      const double y1, const double y2)
    {
        double a, b, c, d;
        double ans;
        int sgn = sign(xy0);
   
        int config = sx0+2*sx1+4*sy0+8*sy1;
        switch(config) {
        case 0:   
        case 2:
        case 8:
        case 10:
            std::cerr << "Oops! case " << d << " in TValFormula\n\n";
            exit(1);
        case 1:
        case 9:
//      sgn = x1==0. ? sgn : sign(x1);
            a = 1./dx/dx;
            b = -2.*x1/dx/dx;
            c = -1.+sqr(x1/dx);
            ans = (-b+sgn*sqrt(b*b-4*a*c))/2/a;
            break;
        case 3:
        case 11:
//      sgn = x1==0. ? sgn : sign(x1);
            a = 2.25/dx/dx;
            b = (-6*x1+1.5*x2)/dx/dx;
            c = -1+sqr((2*x1-x2/2)/dx);
            ans = (-b+sgn*sqrt(b*b-4*a*c))/2/a;
            break;
        case 4:
        case 6:
//      sgn = y1==0. ? sgn : sign(y1);
            a = 1./dy/dy;
            b = -2.*y1/dy/dy;
            c = -1+sqr(y1/dy);
            ans = (-b+sgn*sqrt(b*b-4*a*c))/2/a;
            break;
        case 5:
//      sgn = (x1+y1==0.) ? sgn : sign(x1+y1);
            a = sqr(1/dx)+sqr(1/dy);
            b = -2*(x1/dx/dx+y1/dy/dy);
            c = -1+sqr(x1/dx)+sqr(y1/dy);
            d = b*b-4*a*c;
            if (d < 0) {
                if (fabs(x1) < fabs(y1)) {
                    ans = TValFormula(1,0,0,0,xy0,x1,x2,y1,y2);
                }
                else {
                    ans = TValFormula(0,0,1,0,xy0,x1,x2,y1,y2);
                }
            }
            else
                ans = (-b+sgn*sqrt(d))/2/a;
            break;
        case 7:
//      sgn = (x1+y1==0.) ? sgn : sign(x1+y1);
            a = 2.25*sqr(1/dx)+sqr(2/dy);
            b = (-6*x1+1.5*x2)/dx/dx-2*y1/dy/dy;
            c = -1+sqr((2*x1-x2/2)/dx)+sqr(y1/dy);
            d = b*b-4*a*c;
            if (d < 0) {
                if (fabs(x1) < fabs(y1)) {
                    ans = TValFormula(1,1,0,0,xy0,x1,x2,y1,y2);
                }
                else {
                    ans = TValFormula(0,0,1,0,xy0,x1,x2,y1,y2);
                }
            }
            else
                ans = (-b+sgn*sqrt(d))/2/a;
            break;
        case 12:
        case 14:
//      sgn = y1==0. ? sgn : sign(y1);
            a = 2.25/dy/dy;
            b = (-6*y1+1.5*y2)/dy/dy;
            c = -1+sqr((2*y1-y2/2)/dy);
            ans = (-b+sgn*sqrt(b*b-4*a*c))/2/a;
            break;
        case 13:
//      sgn = (x1+y1==0.) ? sgn : sign(x1+y1);
            a = sqr(1/dx)+sqr(1.5/dy);
            b = -2*x1/dx/dx+(-6*y1+1.5*y2)/dy/dy/2.;
            c = -1+sqr(x1/dx)+sqr((2*y1-y2/2)/dy);
            d = b*b-4*a*c;
            if (d < 0) {
                if (fabs(x1) < fabs(y1)) {
                    ans = TValFormula(1,0,0,0,xy0,x1,x2,y1,y2);
                }
                else {
                    ans = TValFormula(0,0,1,1,xy0,x1,x2,y1,y2);
                }
            }
            else
                ans = (-b+sgn*sqrt(d))/2/a;
            break;
        case 15:
//      sgn = (x1+y1==0.) ? sgn : sign(x1+y1);
            a = sqr(1.5/dx)+sqr(1.5/dy);
            b = (-6*x1+1.5*x2)/dx/dx+(-6*y1+1.5*y2)/dy/dy;
            c = -1+sqr((2*x1-x2/2)/dx)+sqr((2*y1-y2/2)/dy);
            d = b*b-4*a*c;
            if (d < 0) {
                if (fabs(x1) < fabs(y1)) {
                    ans = TValFormula(1,1,0,0,xy0,x1,x2,y1,y2);
                }
                else {
                    ans = TValFormula(0,0,1,1,xy0,x1,x2,y1,y2);
                }
            }
            else
                ans = (-b+sgn*sqrt(d))/2/a;
            break;
        }
   
        return ans;
    }
   
    double UniformMesh2D::ComputeTVal(const int i, const int j, const int k,
                                      const int ikus)
    {
        int sx[2][2];
        int sy[2][2];

        sx[0][0] = idata_(i-1,j,ikus) == 2 ? 1 : 0;
        sx[0][1] = (idata_(i-2,j,ikus) == 2 && sx[0][0]) ? 1 : 0;
        sx[1][0] = idata_(i+1,j,ikus) == 2 ? 1 : 0;
        sx[1][1] = (idata_(i+2,j,ikus) == 2 && sx[1][0]) ? 1 : 0;
        sy[0][0] = idata_(i,j-1,ikus) == 2 ? 1 : 0;
        sy[0][1] = (idata_(i,j-2,ikus) == 2 && sy[0][0]) ? 1 : 0;
        sy[1][0] = idata_(i,j+1,ikus) == 2 ? 1 : 0;
        sy[1][1] = (idata_(i,j+2,ikus) == 2 && sy[1][0]) ? 1 : 0;
   
        double temptval;
        double tval = DBL_MAX;

        // Down left

        if (sx[0][0] || sy[0][0]) {
            temptval = TValFormula(sx[0][0],sx[0][1],sy[0][0],sy[0][1],
                                   data_(i,j,k),data_(i-1,j,k),data_(i-2,j,k),
                                   data_(i,j-1,k),data_(i,j-2,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
   
        // Down right

        if (sx[1][0] || sy[0][0]) {
            temptval = TValFormula(sx[1][0],sx[1][1],sy[0][0],sy[0][1],
                                   data_(i,j,k),data_(i+1,j,k),data_(i+2,j,k),
                                   data_(i,j-1,k),data_(i,j-2,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }

// Up left
   
        if (sx[0][0] || sy[1][0]) {
            temptval = TValFormula(sx[0][0],sx[0][1],sy[1][0],sy[1][1],
                                   data_(i,j,k),data_(i-1,j,k),data_(i-2,j,k),
                                   data_(i,j+1,k),data_(i,j+2,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
   
// Up right
   
        if (sx[1][0] || sy[1][0]) {
            temptval = TValFormula(sx[1][0],sx[1][1],sy[1][0],sy[1][1],
                                   data_(i,j,k),data_(i+1,j,k),data_(i+2,j,k),
                                   data_(i,j+1,k),data_(i,j+2,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
   
        return data_(i,j,k) = tval;
    }
    
	 //new periodic versions of all functions for the velocity extensions
    inline int mod(const int a, const int b) {return (a%b+b)%b;}
	 
	 void UniformMesh2D::ExtendVelocity(const int ktemp, const int kvel,
                                       const int ikus, const int ikhi,
                                       const int imask, const bool xper, const bool yper)
    {  
        NumericTrait<double> t;
      
        // Initialize the accepted and tentative points
      
        // Initialize the UPDATE_STATUS

        int i, j;
        int iP, iM, jP, jM;
		  for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                    data_(i,j,ktemp) = 0.;
                } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,ktemp) = t.max_value;
                }
   
        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);
		  
		  
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 2) {
                    /* original
						  if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                    if (i < maxi-1)
                        UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                    if (j < maxj-1)
                        UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                    */
						  iP = xper ? mod(i+1,maxi-1) : i+1;
						  jP = yper ? mod(j+1,maxj-1) : j+1;
						  iM = xper ? mod(i-1,maxi-1) : i-1;
						  jM = yper ? mod(j-1,maxj-1) : j-1;
						  if (iM >= 0) UpdateNeighbor(theHeap, iM, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                    if (jM >= 0) UpdateNeighbor(theHeap, i, jM, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                    if (iP <= maxi-1)
                        UpdateNeighbor(theHeap, iP, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                    if (jP <= maxj-1)
                        UpdateNeighbor(theHeap, i, jP, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                
					 
					 }
            
        {
            int i, j;
				int iP, jP, iM, jM;
//         PlotWindow2D exdisp(true);
//         exdisp.SetRange(-1., maxi, -1., maxj);
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
//            exdisp.Dot(i,j);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask))
                    ExtendVals(i, j, ktemp, ikus, ikhi, &kvel, 1, xper, yper);

					  iP = xper ? mod(i+1,maxi-1) : i+1;
					  jP = yper ? mod(j+1,maxj-1) : j+1;
					  iM = xper ? mod(i-1,maxi-1) : i-1;
					  jM = yper ? mod(j-1,maxj-1) : j-1;
					 
					 /* original
					 if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
            	 */
					  if (iM >= 0) UpdateNeighbor(theHeap, iM, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
					  if (jM >= 0) UpdateNeighbor(theHeap, i, jM, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
					  if (iP <= maxi-1)
							UpdateNeighbor(theHeap, iP, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
					  if (jP <= maxj-1)
							UpdateNeighbor(theHeap, i, jP, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
				}
        }

    }      
	 //overloaded for initial velocities at non-zero phi values
	 void UniformMesh2D::ExtendVelocity(const int kphi, const int ktemp, const int kvel,
                                       const int ikus, const int ikhi,
                                       const int imask, const bool xper, const bool yper)
    {  
        NumericTrait<double> t;
      
        // Initialize the accepted and tentative points
      
        // Initialize the UPDATE_STATUS

        int i, j;
        int iP, iM, jP, jM;
		  for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                    //data_(i,j,ktemp) = 0.;
                    data_(i,j,ktemp) = fabs(data_(i,j,kphi));
					 } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,ktemp) = t.max_value;
                	  //data_(i,j,ktemp) = fabs(data_(i,j,kphi));
					 }
		  
		  
		  WriteBinary("phiCpyBefore.out", ktemp);

        // Now do fast marching to complete extension

        UniformMesh2D_Heap theHeap(this);
		  
		  
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 2) {
                    /* original
						  if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                    if (i < maxi-1)
                        UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                    if (j < maxj-1)
                        UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                    */
						  iP = xper ? mod(i+1,maxi-1) : i+1;
						  jP = yper ? mod(j+1,maxj-1) : j+1;
						  iM = xper ? mod(i-1,maxi-1) : i-1;
						  jM = yper ? mod(j-1,maxj-1) : j-1;
						  if (iM >= 0) UpdateNeighbor(theHeap, iM, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                    if (jM >= 0) UpdateNeighbor(theHeap, i, jM, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                    if (iP <= maxi-1)
                        UpdateNeighbor(theHeap, iP, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                    if (jP <= maxj-1)
                        UpdateNeighbor(theHeap, i, jP, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
                
					 
					 }
            
        {
            int i, j;
				int iP, jP, iM, jM;
//         PlotWindow2D exdisp(true);
//         exdisp.SetRange(-1., maxi, -1., maxj);
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
//            exdisp.Dot(i,j);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask))
                    ExtendVals(i, j, ktemp, ikus, ikhi, &kvel, 1, xper, yper);

					  iP = xper ? mod(i+1,maxi-1) : i+1;
					  jP = yper ? mod(j+1,maxj-1) : j+1;
					  iM = xper ? mod(i-1,maxi-1) : i-1;
					  jM = yper ? mod(j-1,maxj-1) : j-1;
					 
					 /* original
					 if (i > 0) UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
            	 */
					  if (iM >= 0) UpdateNeighbor(theHeap, iM, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
					  if (jM >= 0) UpdateNeighbor(theHeap, i, jM, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
					  if (iP <= maxi-1)
							UpdateNeighbor(theHeap, iP, j, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
					  if (jP <= maxj-1)
							UpdateNeighbor(theHeap, i, jP, ktemp, ikus, ikhi, xper, yper, data_(i,j,ktemp));
				}
        }

    }      
    
	 void UniformMesh2D::UpdateNeighbor(UniformMesh2D_Heap& heap, const int i,
                                       const int j, const int k, const int ikus,
                                       const int ikhi, const bool xper, const bool yper,
													const double phiAcpt)
    {
#ifdef LEVEL_DEBUG
        int ii;
        double d;
#endif
        switch(idata_(i,j,ikus)) {
        case 0:
            heap.Insert(i,j,fabs(ComputeTVal(i,j,k,ikus,xper,yper,phiAcpt)),ikhi);
            idata_(i,j,ikus) = 1;
            break;
        case 1:
#ifdef LEVEL_DEBUG
            ii = idata_(i,j,ikhi);
            d = fabs(ComputeTVal(i,j,k,ikus,xper,yper,phiAcpt));
            heap.Change(ii,d,ikhi);
#else
            heap.Change(idata_(i,j,ikhi),fabs(ComputeTVal(i,j,k,ikus,xper,yper,phiAcpt)),ikhi);
#endif
            break;
        case 2:
            break;
        default:
            std::cerr << "Error: Gridpoint " << i << "," << j << " not initialized.\n";
            break;
        }
    }
    
    double UniformMesh2D::ComputeTVal(const int i, const int j, const int k,
                                      const int ikus, const bool xper, const bool yper, const double phiAcpt)
    {
        int sx[2][2];
        int sy[2][2];
		  
		  int iP, jP, iM, jM, iPP, jPP, iMM, jMM;
		  
		  iP = xper ? mod(i+1,maxi-1) : i+1;
		  jP = yper ? mod(j+1,maxj-1) : j+1;
		  iM = xper ? mod(i-1,maxi-1) : i-1;
		  jM = yper ? mod(j-1,maxj-1) : j-1;
		  iPP = xper ? mod(i+2,maxi-1) : i+2;
		  jPP = yper ? mod(j+2,maxj-1) : j+2;
		  iMM = xper ? mod(i-2,maxi-1) : i-2;
		  jMM = yper ? mod(j-2,maxj-1) : j-2;
		  
		  /* original
        sx[0][0] = idata_(i-1,j,ikus) == 2 ? 1 : 0;
        sx[0][1] = (idata_(i-2,j,ikus) == 2 && sx[0][0]) ? 1 : 0;
        sx[1][0] = idata_(i+1,j,ikus) == 2 ? 1 : 0;
        sx[1][1] = (idata_(i+2,j,ikus) == 2 && sx[1][0]) ? 1 : 0;
        sy[0][0] = idata_(i,j-1,ikus) == 2 ? 1 : 0;
        sy[0][1] = (idata_(i,j-2,ikus) == 2 && sy[0][0]) ? 1 : 0;
        sy[1][0] = idata_(i,j+1,ikus) == 2 ? 1 : 0;
        sy[1][1] = (idata_(i,j+2,ikus) == 2 && sy[1][0]) ? 1 : 0;
        */
        
		  sx[0][0] = (idata_(iM,j,ikus) == 2 && data_(iM,j,k) <= phiAcpt ) ? 1 : 0;
        sx[0][1] = ((idata_(iMM,j,ikus) == 2 && sx[0][0]) && data_(iMM,j,k) <= phiAcpt ) ? 1 : 0;
        sx[1][0] = (idata_(iP,j,ikus) == 2 && data_(iP,j,k) <= phiAcpt) ? 1 : 0;
        sx[1][1] = ((idata_(iPP,j,ikus) == 2 && sx[1][0]) && data_(iPP,j,k) <= phiAcpt ) ? 1 : 0;
        sy[0][0] = (idata_(i,jM,ikus) == 2 && data_(i,jM,k) <= phiAcpt ) ? 1 : 0;
        sy[0][1] = ((idata_(i,jMM,ikus) == 2 && sy[0][0]) && data_(i,jMM,k) <= phiAcpt ) ? 1 : 0;
        sy[1][0] = (idata_(i,jP,ikus) == 2 && data_(i,jP,k) <= phiAcpt ) ? 1 : 0;
        sy[1][1] = ((idata_(i,jPP,ikus) == 2 && sy[1][0]) && data_(i,jPP,k) <= phiAcpt ) ? 1 : 0;
		 	
		 //DEBUG (try forcing it to do first order only?)
		 sx[0][1] = 0;
		 sx[1][1] = 0;
		 sy[0][1] = 0;
		 sy[1][1] = 0;
		  
		  double temptval;
        double tval = DBL_MAX;

        // Down left

        if (sx[0][0] || sy[0][0]) {
            /*
				temptval = TValFormula(sx[0][0],sx[0][1],sy[0][0],sy[0][1],
                                   data_(i,j,k),data_(i-1,j,k),data_(i-2,j,k),
                                   data_(i,j-1,k),data_(i,j-2,k));
            */
				temptval = TValFormula(sx[0][0],sx[0][1],sy[0][0],sy[0][1],
                                   data_(i,j,k),data_(iM,j,k),data_(iMM,j,k),
                                   data_(i,jM,k),data_(i,jMM,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
   
        // Down right

        if (sx[1][0] || sy[0][0]) {
            /*
				temptval = TValFormula(sx[1][0],sx[1][1],sy[0][0],sy[0][1],
                                   data_(i,j,k),data_(i+1,j,k),data_(i+2,j,k),
                                   data_(i,j-1,k),data_(i,j-2,k));
            */
				temptval = TValFormula(sx[1][0],sx[1][1],sy[0][0],sy[0][1],
                                   data_(i,j,k),data_(iP,j,k),data_(iPP,j,k),
                                   data_(i,jM,k),data_(i,jMM,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }

// Up left
   
        if (sx[0][0] || sy[1][0]) {
            /*
				temptval = TValFormula(sx[0][0],sx[0][1],sy[1][0],sy[1][1],
                                   data_(i,j,k),data_(i-1,j,k),data_(i-2,j,k),
                                   data_(i,j+1,k),data_(i,j+2,k));
            */
				temptval = TValFormula(sx[0][0],sx[0][1],sy[1][0],sy[1][1],
                                   data_(i,j,k),data_(iM,j,k),data_(iMM,j,k),
                                   data_(i,jP,k),data_(i,jPP,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
   
// Up right
   
        if (sx[1][0] || sy[1][0]) {
            /*
				temptval = TValFormula(sx[1][0],sx[1][1],sy[1][0],sy[1][1],
                                   data_(i,j,k),data_(i+1,j,k),data_(i+2,j,k),
                                   data_(i,j+1,k),data_(i,j+2,k));
            */
				temptval = TValFormula(sx[1][0],sx[1][1],sy[1][0],sy[1][1],
                                   data_(i,j,k),data_(iP,j,k),data_(iPP,j,k),
                                   data_(i,jP,k),data_(i,jPP,k));
            tval = fabs(temptval) < fabs(tval) ? temptval : tval;
        }
   
        return data_(i,j,k) = tval;
    }
    
	 void UniformMesh2D::ExtendVals(const int i, const int j, const int k,
                                   const int ikus, const int ikhi,
                                   const int* p, const int pl, const bool xper, const bool yper)
    {
        
		  int iP, jP, iM, jM;
		  iP = xper ? mod(i+1,maxi-1) : i+1;
		  jP = yper ? mod(j+1,maxj-1) : j+1;
		  iM = xper ? mod(i-1,maxi-1) : i-1;
		  jM = yper ? mod(j-1,maxj-1) : j-1;
		  
		  
		  int d;
        //int sgn = data_(i,j,k) >= 0. ? 1 : -1;
        /* original
		  d  = ((i > 0) && idata_(i-1,j,ikus)==2) ? 1 : 0;
        d += ((j > 0) && idata_(i,j-1,ikus)==2) ? 2 : 0;
        d += ((i < maxi-1) && idata_(i+1,j,ikus)==2) ? 4 : 0;
        d += ((j < maxj-1) && idata_(i,j+1,ikus)==2) ? 8 : 0;
        */
		  
		  d  = ((iM >= 0) && idata_(iM,j,ikus)==2) ? 1 : 0;
        d += ((jM >= 0) && idata_(i,jM,ikus)==2) ? 2 : 0;
        d += ((iP <= maxi-1) && idata_(iP,j,ikus)==2) ? 4 : 0;
        d += ((jP <= maxj-1) && idata_(i,jP,ikus)==2) ? 8 : 0;
   
        double a;
        double b;
        int l;
		  /* original
        switch(d) {
        case 0:
            std::cout << "Case 0 hit\n";
            a = (data_(i+1,j,k)+data_(i,j+1,k)+data_(i-1,j,k)+data_(i,j-1,k)
                 -4.*data_(i,j,k))/4.;
            a = 0.;
            b = 0.;
            for (l=0; l<pl; ++l) {
                if (idata_(i+1,j,ikhi)) {
                    a += data_(i+1,j,p[l]);
                    b += 1.;
                }
                if (idata_(i-1,j,ikhi)) {
                    a += data_(i-1,j,p[l]);
                    b += 1.;
                }
                if (idata_(i,j+1,ikhi)) {
                    a += data_(i,j+1,p[l]);
                    b += 1.;
                }
                if (idata_(i,j-1,ikhi)) {
                    a += data_(i,j-1,p[l]);
                    b += 1.;
                }
                data_(i,j,p[l]) = a/b;
            }
            break;
        case 1:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i-1,j,p[l]);
            break;
        case 2:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i,j-1,p[l]);
            break;
        case 3:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)
                                          -data_(i,j-1,k))*data_(i,j-1,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(i-1,j,k))
                                   *data_(i-1,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                      +dy*dy*(data_(i,j,k)-data_(i-1,j,k)));
            break;
        case 4:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i+1,j,p[l]);
            break;
        case 5:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i-1,j,p[l])*(data_(i,j,k)-data_(i-1,j,k))
                                   -data_(i+1,j,p[l])*(data_(i+1,j,k)-data_(i,j,k)))
                    /(2.*data_(i,j,k)-data_(i-1,j,k)-data_(i+1,j,k));
            break;
        case 6:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                                   *data_(i,j-1,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(i+1,j,k))
                                   *data_(i+1,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                      +dy*dy*(data_(i,j,k)-data_(i+1,j,k)));
            break;
        case 7:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                                   *data_(i,j-1,p[l])
                                   +dy*dy*(data_(i-1,j,p[l])
                                           *(data_(i,j,k)-data_(i-1,j,k))
                                           -data_(i+1,j,p[l])
                                           *(data_(i+1,j,k)-data_(i,j,k))))
                    /(dx*dx*(data_(i,j,k)-data_(i,j-1,k))
                      +dy*dy*(2*data_(i,j,k)-data_(i-1,j,k)
                              -data_(i+1,j,k)));
            break;
        case 8:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i,j+1,p[l]);
            break;
        case 9:
            for (l=0; l<pl; ++l) {
                a = dx*dx*(data_(i,j,k)-data_(i,j+1,k));
                b = dy*dy*(data_(i,j,k)-data_(i-1,j,k));
                data_(i,j,p[l]) = (a*data_(i,j+1,p[l])+b*data_(i-1,j,p[l]))/(a+b);
            }
            break;
        case 10:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i,j-1,p[l])*(data_(i,j,k)-data_(i,j-1,k))
                                   -data_(i,j+1,p[l])*(data_(i,j+1,k)
                                                       -data_(i,j,k)))
                    /(2*data_(i,j,k)-data_(i,j-1,k)-data_(i,j+1,k));
            break;
        case 11:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j+1,p[l])
                                          *(data_(i,j,k)-data_(i,j+1,k))
                                          +data_(i,j-1,p[l])
                                          *(data_(i,j,k)-data_(i,j-1,k)))
                                   +dy*dy*data_(i-1,j,p[l])
                                   *(data_(i,j,k)-data_(i-1,j,k)))
                    /(dx*dx*(2*data_(i,j,k)
                             -data_(i,j-1,k)-data_(i,j+1,k))
                      +dy*dy*(data_(i,j,k)-data_(i-1,j,k)));
            break;
        case 12:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)
                                          -data_(i,j+1,k))*data_(i,j+1,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(i+1,j,k))
                                   *data_(i+1,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,j+1,k))
                      +dy*dy*(data_(i,j,k)-data_(i+1,j,k)));
            break;
        case 13:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*data_(i,j+1,p[l])
                                   *(data_(i,j,k)-data_(i,j+1,k))
                                   +dy*dy*(data_(i-1,j,p[l])
                                           *(data_(i,j,k)-data_(i-1,j,k))
                                           +data_(i+1,j,p[l])
                                           *(data_(i,j,k)-data_(i+1,j,k))))
                    /(dx*dx*(data_(i,j,k)-data_(i,j+1,k))
                      +dy*dy*(2*data_(i,j,k)
                              -data_(i-1,j,k)-data_(i+1,j,k)));
            break;
        case 14:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j-1,p[l])
                                          *(data_(i,j,k)-data_(i,j-1,k))
                                          +data_(i,j+1,p[l])
                                          *(data_(i,j,k)-data_(i,j+1,k)))
                                   +dy*dy*data_(i+1,j,p[l])
                                   *(data_(i,j,k)-data_(i+1,j,k)))
                    /(dx*dx*(2*data_(i,j,k)-data_(i,j-1,k)
                             -data_(i,j+1,k))
                      +dy*dy*(data_(i,j,k)-data_(i+1,j,k)));
            break;
        case 15:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i-1,j,p[l])+data_(i,j-1,p[l])
                                   +data_(i+1,j,p[l])+data_(i,j+1,p[l]))/4.;
            break;
        }
		  */
        
		  switch(d) {
        case 0:
            std::cout << "Case 0 hit\n";
            a = (data_(iP,j,k)+data_(i,jP,k)+data_(iM,j,k)+data_(i,jM,k)
                 -4.*data_(i,j,k))/4.;
            a = 0.;
            b = 0.;
            for (l=0; l<pl; ++l) {
                if (idata_(iP,j,ikhi)) {
                    a += data_(iP,j,p[l]);
                    b += 1.;
                }
                if (idata_(iM,j,ikhi)) {
                    a += data_(iM,j,p[l]);
                    b += 1.;
                }
                if (idata_(i,jP,ikhi)) {
                    a += data_(i,jP,p[l]);
                    b += 1.;
                }
                if (idata_(i,jM,ikhi)) {
                    a += data_(i,jM,p[l]);
                    b += 1.;
                }
                data_(i,j,p[l]) = a/b;
            }
            break;
        case 1:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(iM,j,p[l]);
            break;
        case 2:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i,jM,p[l]);
            break;
        case 3:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)
                                          -data_(i,jM,k))*data_(i,jM,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(iM,j,k))
                                   *data_(iM,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,jM,k))
                      +dy*dy*(data_(i,j,k)-data_(iM,j,k)));
            break;
        case 4:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(iP,j,p[l]);
            break;
        case 5:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(iM,j,p[l])*(data_(i,j,k)-data_(iM,j,k))
                                   -data_(iP,j,p[l])*(data_(iP,j,k)-data_(i,j,k)))
                    /(2.*data_(i,j,k)-data_(iM,j,k)-data_(iP,j,k));
            break;
        case 6:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)-data_(i,jM,k))
                                   *data_(i,jM,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(iP,j,k))
                                   *data_(iP,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,jM,k))
                      +dy*dy*(data_(i,j,k)-data_(iP,j,k)));
            break;
        case 7:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)-data_(i,jM,k))
                                   *data_(i,jM,p[l])
                                   +dy*dy*(data_(iM,j,p[l])
                                           *(data_(i,j,k)-data_(iM,j,k))
                                           -data_(iP,j,p[l])
                                           *(data_(iP,j,k)-data_(i,j,k))))
                    /(dx*dx*(data_(i,j,k)-data_(i,jM,k))
                      +dy*dy*(2*data_(i,j,k)-data_(iM,j,k)
                              -data_(iP,j,k)));
            break;
        case 8:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = data_(i,jP,p[l]);
            break;
        case 9:
            for (l=0; l<pl; ++l) {
                a = dx*dx*(data_(i,j,k)-data_(i,jP,k));
                b = dy*dy*(data_(i,j,k)-data_(iM,j,k));
                data_(i,j,p[l]) = (a*data_(i,jP,p[l])+b*data_(iM,j,p[l]))/(a+b);
            }
            break;
        case 10:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i,jM,p[l])*(data_(i,j,k)-data_(i,jM,k))
                                   -data_(i,jP,p[l])*(data_(i,jP,k)
                                                       -data_(i,j,k)))
                    /(2*data_(i,j,k)-data_(i,jM,k)-data_(i,jP,k));
            break;
        case 11:
            for (l=0; l<pl; ++l) 
					 data_(i,j,p[l]) = (dx*dx*(data_(i,jP,p[l])
                                          *(data_(i,j,k)-data_(i,jP,k))
                                          +data_(i,jM,p[l])
                                          *(data_(i,j,k)-data_(i,jM,k)))
                                   +dy*dy*data_(iM,j,p[l])
                                   *(data_(i,j,k)-data_(iM,j,k)))
                    /(dx*dx*(2*data_(i,j,k)
                             -data_(i,jM,k)-data_(i,jP,k))
                      +dy*dy*(data_(i,j,k)-data_(iM,j,k)));
            break;
        case 12:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,j,k)
                                          -data_(i,jP,k))*data_(i,jP,p[l])
                                   +dy*dy*(data_(i,j,k)-data_(iP,j,k))
                                   *data_(iP,j,p[l]))
                    /(dx*dx*(data_(i,j,k)-data_(i,jP,k))
                      +dy*dy*(data_(i,j,k)-data_(iP,j,k)));
            break;
        case 13:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*data_(i,jP,p[l])
                                   *(data_(i,j,k)-data_(i,jP,k))
                                   +dy*dy*(data_(iM,j,p[l])
                                           *(data_(i,j,k)-data_(iM,j,k))
                                           +data_(iP,j,p[l])
                                           *(data_(i,j,k)-data_(iP,j,k))))
                    /(dx*dx*(data_(i,j,k)-data_(i,jP,k))
                      +dy*dy*(2*data_(i,j,k)
                              -data_(iM,j,k)-data_(iP,j,k)));
            break;
        case 14:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (dx*dx*(data_(i,jM,p[l])
                                          *(data_(i,j,k)-data_(i,jM,k))
                                          +data_(i,jP,p[l])
                                          *(data_(i,j,k)-data_(i,jP,k)))
                                   +dy*dy*data_(iP,j,p[l])
                                   *(data_(i,j,k)-data_(iP,j,k)))
                    /(dx*dx*(2*data_(i,j,k)-data_(i,jM,k)
                             -data_(i,jP,k))
                      +dy*dy*(data_(i,j,k)-data_(iP,j,k)));
            break;
        case 15:
            for (l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(iM,j,p[l])+data_(i,jM,p[l])
                                   +data_(iP,j,p[l])+data_(i,jP,p[l]))/4.;
            break;
        }


    }

	 

} //end of the namespace 

