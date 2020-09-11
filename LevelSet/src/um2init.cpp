#include "defs.h"
#include <math.h>
#include "uniformmesh2d.h"
#include "um2boundary.h"
#ifdef LEVEL_DEBUG
#include <fstream>
#endif

namespace levelset {
        
    void UniformMesh2D::InitInterface(const InitialFunc& f, const double start,
                                      const double end, const int func, const int ktemp,
                                      const int ik1, const int ik2)
    {
        double thresh = sqr(min(dx,dy)/10.);
        double ds = thresh;
        int i[2], j[2];
        double ifrac[2], jfrac[2];
        double x[2], y[2];
        x[0] = f.X(start);
        y[0] = f.Y(start);
        LocToIndex(x[0], y[0], i[0], j[0], ifrac[0], jfrac[0]);
        for (int k=0; k<3; ++k)
            ds = thresh*ds/sqrt(sqr(f.X(start+ds)-x[0])+sqr(f.Y(start+ds)-y[0]));
   
        for (double s=start; s<=end; s+=ds) {
            ds = thresh*ds/sqrt(sqr(f.X(s+ds)-x[0])+sqr(f.Y(s+ds)-y[0]));
            x[1] = f.X(s+ds);
            y[1] = f.Y(s+ds);
            LocToIndex(x[1], y[1], i[1], j[1], ifrac[1], jfrac[1]);
            data_(i[0], j[0], func) = sign(data_(i[0],j[0],func))
                *min(double(sqrt(sqr(ifrac[0]*dx)+sqr(jfrac[0]*dy))),
                     fabs(data_(i[0],j[0],func)));
            data_(i[0]+1, j[0], func) = sign(data_(i[0]+1,j[0],func))
                *min(double(sqrt(sqr((1-ifrac[0])*dx)+sqr(jfrac[0]*dy))),
                     fabs(data_(i[0]+1,j[0],func)));
            data_(i[0], j[0]+1, func) = sign(data_(i[0],j[0]+1,func))
                *min(double(sqrt(sqr(ifrac[0]*dx)+sqr((1-jfrac[0])*dy))),
                     fabs(data_(i[0],j[0]+1,func)));
            data_(i[0]+1, j[0]+1, func) = sign(data_(i[0]+1,j[0]+1,func))
                *min(double(sqrt(sqr((1-ifrac[0])*dx)+sqr((1-jfrac[0])*dy))),
                     fabs(data_(i[0]+1,j[0]+1,func)));
            if (i[0] != i[1] || (x[0]-zero[0])*(x[1]-zero[0]) < 0.) {
                int ii, jj;
                for (jj=j[0]+1, ii=max(i[0],i[1]); jj<maxj; ++jj)
                    data_(ii,jj,func) *= -1;
            }
            if ((y[0]-zero[1])*(y[1]-zero[1]) < 0.) {
                int ii, jj;
                for (ii=max(i[0],i[1]); ii<maxi; ++ii)
                    for (jj=0; jj<maxj; ++jj)
                        data_(ii,jj,func) *= -1;
            }
            i[0] = i[1];
            j[0] = j[1];
            x[0] = x[1];
            y[0] = y[1];
            ifrac[0] = ifrac[1];
            jfrac[0] = jfrac[1];
        }
        Reinitialize(func, ktemp, ik1, ik2);
    }
   
    void UniformMesh2D::SetValues(const InitialFunc& f, const int wi, const double mult)
    {
        double x;
        double y;

        for (int i=0; i<maxi; ++i) {
            x = i*dx+zero[0];
            for (int j=0; j<maxj; ++j) {
                y = j*dy+zero[1];
                data_(i,j,wi) = f.XY(x,y)*mult;
            }
        }
        bc->Apply(wi);
    }


    void UniformMesh2D::InitInterface(double (*fx)(double), double (*fy)(double),
                                      const double start, const double end,
                                      const int func)
    {
        double thresh = sqr(min(dx,dy)/100.);
        double ds = (end-start)/2./maxi;
        int i[2], j[2];
        double ifrac[2], jfrac[2];
        double x[2], y[2];
        x[0] = (*fx)(start);
        y[0] = (*fy)(start);
        LocToIndex(x[0], y[0], i[0], j[0], ifrac[0], jfrac[0]);
        do {
            ds /= 2.;
            x[1] = (*fx)(start+ds);
            y[1] = (*fy)(start+ds);
            LocToIndex(x[1], y[1], i[1], j[1], ifrac[1], jfrac[1]);
        } while (sqr(x[0]-x[1])+sqr(y[0]-y[1]) > thresh);
   
        for (double s=start; s<=end; s+=ds) {
            data_(i[0], j[0], func) = sign(data_(i[0],j[0],func))
                *min(sqrt(sqr(ifrac[0]*dx)+sqr(jfrac[0]*dy)),
                     fabs(data_(i[0],j[0],func)));
            data_(i[0]+1, j[0], func) = sign(data_(i[0]+1,j[0],func))
                *min(sqrt(sqr((1-ifrac[0])*dx)+sqr(jfrac[0]*dy)),
                     fabs(data_(i[0]+1,j[0],func)));
            data_(i[0], j[0]+1, func) = sign(data_(i[0],j[0]+1,func))
                *min(sqrt(sqr(ifrac[0]*dx)+sqr((1-jfrac[0])*dy)),
                     fabs(data_(i[0],j[0]+1,func)));
            data_(i[0]+1, j[0]+1, func) = sign(data_(i[0]+1,j[0]+1,func))
                *min(sqrt(sqr((1-ifrac[0])*dx)+sqr((1-jfrac[0])*dy)),
                     fabs(data_(i[0]+1,j[0]+1,func)));
            if (i[0] != i[1] || (x[0]-zero[0])*(x[1]-zero[0]) < 0.) {
                int ii, jj;
                for (jj=j[0]+1, ii=max(i[0],i[1]); jj<maxj; ++jj)
                    data_(ii,jj,func) *= -1;
            }
            if ((y[0]-zero[1])*(y[1]-zero[1]) < 0.) {
                int ii, jj;
                for (ii=max(i[0],i[1]); ii<maxi; ++ii)
                    for (jj=0; jj<maxj; ++jj)
                        data_(ii,jj,func) *= -1;
            }
            x[0] = x[1];
            y[0] = y[1];
            ds *= 4.;
            do {
                ds /= 2.;
                x[1] = (*fx)(s+ds);
                y[1] = (*fy)(s+ds);
                LocToIndex(x[1], y[1], i[1], j[1], ifrac[1], jfrac[1]);
            } while (sqr(x[0]-x[1])+sqr(y[0]-y[1]) > thresh);
        }
        bc->Apply(func);
    }
   
    void UniformMesh2D::SetValues(double (*f)(double,double), const int wi, const double mult)
    {
        double x;
        double y;

        for (int i=0; i<maxi; ++i) {
            x = i*dx+zero[0];
            for (int j=0; j<maxj; ++j) {
                y = j*dy+zero[1];
                data_(i,j,wi) = (*f)(x,y)*mult;
            }
        }
        bc->Apply(wi);
    }

    void UniformMesh2D::SetValues(const double* f, const int wi, const double mult)
    {
        for (int i=0; i<maxi; ++i) {
            for (int j=0; j<maxj; ++j) {
                data_(i,j,wi) = f[i*maxj+j]*mult;
            }
        }
        bc->Apply(wi);
    }

    void UniformMesh2D::SetValues(const double c, const int wi)
    {
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j)
                data_(i,j,wi) = c;
        bc->Apply(wi);
    }

}
