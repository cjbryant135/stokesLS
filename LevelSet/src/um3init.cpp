#include "defs.h"
#include <math.h>
#include "uniformmesh3d.h"
#include "um3boundary.h"
#ifdef LEVEL_DEBUG
#include <fstream>
#endif

namespace levelset {
        
    void UniformMesh3D::InitInterface(const InitialFunc3D& f, const int func, 
                                      const int ktemp, const int ik1, const int ik2)
    {
        double x, y, z;
   
        for (int i=0; i<maxi; ++i) {
            x = i*dx+zero[0];
            for (int j=0; j<maxj; ++j) {
                y = j*dy+zero[1];
                for (int k=0; k<maxk; ++k) {
                    z = k*dz+zero[2];
                    data_(i,j,k,func) = f.XYZ(x,y,z);
                }
            }
        }
            
        bc->Apply(func);
        Reinitialize(func, ktemp, ik1, ik2);
    }
   
    void UniformMesh3D::SetValues(const InitialFunc3D& f, const int func)
    {
        double x, y, z;
   
        for (int i=0; i<maxi; ++i) {
            x = X(i);
            for (int j=0; j<maxj; ++j) {
                y = Y(j);
                for (int k=0; k<maxk; ++k) {
                    z = Z(k);
                    data_(i,j,k,func) = f.XYZ(x,y,z);
                }
            }
        }
            
        bc->Apply(func);
    }


    void UniformMesh3D::InitInterface(double (*fx)(double,double), double (*fy)(double,double),
                                      double (*fz)(double,double), const double start1, 
                                      const double end1, const double start2, 
                                      const double end2, const int func)
    {
        long numsteps = long(max((start1-end1)/min(dx,dy,dz)*4.,
                                 (start2-end2)/min(dx,dy,dz)*4.));
        double ds = (end1-start1)/numsteps;
        double dt = (end2-start2)/numsteps;
        double s, t;
        double x[4], y[4], z[4];
        double mx, my, mz;
        double rx, ry, rz;

        for (int i=0; i<maxi; ++i) {
            rx = X(i);
            for (int j=0; j<maxj; ++j) {
                ry = Y(j);
                for (int k=0; k<maxk; ++k) {
                    rz = Z(k);
            
                    data_(i,j,k,func) = DBL_MAX;
            
                    for (int l=0; l<numsteps; ++l) {
                        s = start1 + l*ds;
                        for (int m=0; m<numsteps; ++m) {
                            t = start2 + m*dt;
      
                            x[0] = fx(s, t);
                            y[0] = fy(s, t);
                            z[0] = fz(s, t);
                            x[1] = fx(s+ds, t);
                            y[1] = fy(s+ds, t);
                            z[1] = fz(s+ds, t);
                            x[2] = fx(s, t+dt);
                            y[2] = fy(s, t+dt);
                            z[2] = fz(s, t+dt);
                            x[3] = fx(s+ds, t+dt);
                            y[3] = fy(s+ds, t+dt);
                            z[3] = fz(s+ds, t+dt);

                            mx = (x[0]+x[1]+x[2]+x[3])/4;
                            my = (y[0]+y[1]+y[2]+y[3])/4;
                            mz = (z[0]+z[1]+z[2]+z[3])/4;
         
                            data_(i,j,k,func) = min(min(data_(i,j,k,func), 
                                                        sqrt(sqr(rx-mx)+sqr(ry-my)+sqr(rz-mz)),
                                                        sqrt(sqr(rx-x[0])+sqr(ry-y[0])
                                                             +sqr(rz-z[0]))),
                                                    min(sqrt(sqr(rx-x[1])+sqr(ry-y[1])
                                                             +sqr(rz-z[1])),
                                                        sqrt(sqr(rx-x[2])+sqr(ry-y[2])
                                                             +sqr(rz-z[2])),
                                                        sqrt(sqr(rx-x[3])+sqr(ry-y[3])
                                                             +sqr(rz-z[3]))));
                        }
                    }
                                                   
                }
            }
        }


        for (int l=0; l<numsteps; ++l) {
            s = start1 + l*ds;
            for (int m=0; m<numsteps; ++m) {
                t = start2 + m*dt;
      
                x[0] = fx(s, t);
                y[0] = fy(s, t);
                z[0] = fz(s, t);
                x[1] = fx(s+ds, t);
                y[1] = fy(s+ds, t);
                z[1] = fz(s+ds, t);
                x[2] = fx(s, t+dt);
                y[2] = fy(s, t+dt);
                z[2] = fz(s, t+dt);
                x[3] = fx(s+ds, t+dt);
                y[3] = fy(s+ds, t+dt);
                z[3] = fz(s+ds, t+dt);

                int iimin = min(I(x[0]), I(x[1]), I(x[2]))+1;
                int iimax = max(I(x[0]), I(x[1]), I(x[2]))+1;
                int jjmin = min(J(y[0]), J(y[1]), J(y[2]))+1;
                int jjmax = max(J(y[0]), J(y[1]), J(y[2]))+1;
                double a, b, d;
                for (int i=iimin; i<=iimax; ++i)    
                    for (int j=jjmin; j<=jjmax; ++j) {
                        d = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
                        a = ((X(i)-x[0])*(y[2]-y[0])-(Y(j)-y[0])*(y[1]-y[0]))/d;
                        b = ((Y(j)-y[0])*(x[1]-x[0])-(X(i)-x[0])*(x[2]-x[0]))/d;
                        if (a >= 0 && b >= 0 && a+b <= 1) {
                            for (int k=K(a*(z[1]-z[0])+b*(z[2]-z[0])+z[0])+1; k<maxk; ++k)
                                data_(i,j,k,func) *= -1;
                        }
                    }
            }
        }
   
        bc->Apply(func);
    }
   
    void UniformMesh3D::SetValues(double (*f)(double,double,double),
                                  const int func)
    {
        double x, y, z;

        for (int i=0; i<maxi; ++i) {
            x = X(i);
            for (int j=0; j<maxj; ++j) {
                y = Y(j);
                for (int k=0; k<maxk; ++k) {
                    z = Z(k);
                    data_(i,j,k,func) = (*f)(x,y,z);
                }
            }
        }
        bc->Apply(func);
    }

    void UniformMesh3D::SetValues(const double* f, const int wi)
    {
        for (int i=0; i<maxi; ++i) 
            for (int j=0; j<maxj; ++j)
                for (int k=0; k<maxk; ++k)
                    data_(i,j,k,wi) = f[i*maxj*maxk+j*maxk+k];
      
        bc->Apply(wi);
    }

    void UniformMesh3D::SetValues(const double c, const int wi)
    {
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j)
                for (int k=0; k<maxk; ++k)
                    data_(i,j,k,wi) = c;
        bc->Apply(wi);
    }

}
