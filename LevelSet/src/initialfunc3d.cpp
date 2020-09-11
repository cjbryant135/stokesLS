#include <iostream>
#include <math.h>
#include <string.h>
#include "initialfunc3d.h"
#include <limits.h>
#include <float.h>
#include "defs.h"
#include "utility.h"

namespace levelset {

    void InitialFunc3D::SetParams(const int nparams, const double* params)
    {
        if (nparams > 0) {
            if (parameter) delete[] parameter;
            parameter = new double[nparams];
            for (int i=0; i<nparams; ++i) parameter[i] = params[i];
            param_num = nparams;
        } else {
            if (parameter) delete[] parameter;
            parameter = NULL;
            param_num = 0;
        }
    }

    void InitialFunc3D::SetParams(const double p1)
    {
        if (parameter) delete[] parameter;
        parameter = new double[1];
        parameter[0] = p1;
        param_num = 1;
    }

    void InitialFunc3D::SetParams(const double p1, const double p2)
    {
        if (parameter) delete[] parameter;
        parameter = new double[2];
        parameter[0] = p1;
        parameter[1] = p2;
        param_num = 2;
    }

    void InitialFunc3D::SetParams(const double p1, const double p2,
                                  const double p3)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        parameter[0] = p1;
        parameter[1] = p2;
        parameter[2] = p3;
        param_num = 3;
    }

    void InitialFunc3D::SetParams(const double p1, const double p2,
                                  const double p3, const double p4)
    {
        if (parameter) delete[] parameter;
        parameter = new double[4];
        parameter[0] = p1;
        parameter[1] = p2;
        parameter[2] = p3;
        parameter[3] = p4;
        param_num = 4;
    }

    char InitialFunc3D::in_triangle(const double u1, const double u2,
                                    const double u3, const double v1,
                                    const double v2, const double v3,
                                    const double x1, const double x2,
                                    const double x3) const
    {
        double d = u1*v2-v1*u2;
        if (d != 0.) {
            double w1 = (x1*v2-x2*v1)/d;
            double w2 = (x2*u1-x1*u2)/d;
            if (w1 >= 0. && w2 >= 0. && w1+w2 <= 1) 
                return(x1*(u2*v3-u3*v2)+x2*(u3*v1-u1*v3)+x3*(u1*v2-u2*v1) >= 0. 
                       ? true : false);
            else
                return(false);
        }
        else
            return(false);
    }
            

    double InitialFunc3D::XYZ(const double s, const double t, const double u) const
    {
        double ds = (stop[0]-start[0])/xy_numsteps;
        double dt = (stop[1]-start[1])/xy_numsteps;
        int    cross = 0;
        double val = DBL_MAX;
        double cs, ct;
        double x[4], y[4], z[4];
        double mx, my, mz;

        for (int i=0; i<xy_numsteps; ++i) {
            cs = start[0] + i*ds;
            for (int j=0; j<xy_numsteps; ++j) {
                ct = start[1] + j*dt;
      
                x[0] = X(cs, ct);
                y[0] = Y(cs, ct);
                z[0] = Z(cs, ct);
                x[1] = X(cs+ds, ct);
                y[1] = Y(cs+ds, ct);
                z[1] = Z(cs+ds, ct);
                x[2] = X(cs, ct+dt);
                y[2] = Y(cs, ct+dt);
                z[2] = Z(cs, ct+dt);
                x[3] = X(cs+ds, ct+dt);
                y[3] = Y(cs+ds, ct+dt);
                z[3] = Z(cs+ds, ct+dt);

                mx = (x[0]+x[1]+x[2]+x[3])/4;
                my = (y[0]+y[1]+y[2]+y[3])/4;
                mz = (z[0]+z[1]+z[2]+z[3])/4;
         
                val = int(min(val, sqrt((s-mx)*(s-mx)+(t-my)*(t-my)+(u-mz)*(u-mz))));
         
                if (in_triangle(x[1]-x[0],y[1]-y[0],z[1]-z[0],
                                x[2]-x[0],y[2]-y[0],z[2]-z[0],
                                s-x[0],t-y[0],u-z[0]))
                    ++cross;
                if (in_triangle(x[1]-x[3],y[1]-y[3],z[1]-z[3],
                                x[2]-x[3],y[2]-y[3],z[2]-z[3],
                                s-x[3],t-y[3],u-z[3]))
                    ++cross;
      
            }
        }
        return val*(cross%2 ? 1 : -1);
    }

      
}












