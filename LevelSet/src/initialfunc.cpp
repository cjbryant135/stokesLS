#include <iostream>
#include <math.h>
#include <string.h>
#include "initialfunc.h"
#include <limits.h>
#include <float.h>
#include "defs.h"
#include "utility.h"

namespace levelset {

    void InitialFunc::SetParams(const int nparams, const double* params)
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

    void InitialFunc::SetParams(const double p1)
    {
        if (parameter) delete[] parameter;
        parameter = new double[1];
        parameter[0] = p1;
        param_num = 1;
    }

    void InitialFunc::SetParams(const double p1, const double p2)
    {
        if (parameter) delete[] parameter;
        parameter = new double[2];
        parameter[0] = p1;
        parameter[1] = p2;
        param_num = 2;
    }

    void InitialFunc::SetParams(const double p1, const double p2, const double p3)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        parameter[0] = p1;
        parameter[1] = p2;
        parameter[2] = p3;
        param_num = 3;
    }

    void InitialFunc::SetParams(const double p1, const double p2, const double p3, const double p4)
    {
        if (parameter) delete[] parameter;
        parameter = new double[4];
        parameter[0] = p1;
        parameter[1] = p2;
        parameter[2] = p3;
        parameter[3] = p4;
        param_num = 4;
    }

#if 0
    double InitialFunc::XY(const double s, const double t) const
    {
        double dc = (stop-start)/xy_numsteps;
        double ang = 0;
        double val = DBL_MAX;
        double c;
        double x[2], y[2];
        double mx, my;

        x[0] = X(start);
        y[0] = Y(start);
        for (int i=0; i<xy_numsteps; ++i) {
            c = start+i*dc;
            x[1] = X(c+dc);
            y[1] = Y(c+dc);
            mx = (x[0]+x[1])/2;
            my = (y[0]+y[1])/2;
            val = min(val, sqrt((s-mx)*(s-mx)+(t-my)*(t-my)));
            ang += asin(((y[1]-t)*(x[0]-s)-(y[0]-t)*(x[1]-s))/
                        sqrt((sqr(x[0]-s)+sqr(y[0]-t))
                             *(sqr(x[1]-s)+sqr(y[1]-t))));
            x[0] = x[1];
            y[0] = y[1];
        }
        return val*(abs(ang) > M_PI ? 1 : -1);
    }
#else

    double distfunc(const double s, const double t, const double x, const double y)
    {
        return sqrt((s-x)*(s-x)+(t-y)*(t-y));
    }

    double InitialFunc::Brent(const double bx, const double ds, 
                              const double s, const double t) const
    {
        double a = bx-ds;
        double b = bx+ds;
        double x = bx;
        double w = bx;
        double v = bx;
        double fw = distfunc(s,t,X(x),Y(x));
        double fv = fw;
        double fx = fw;
        double u, fu;
        double e = 0.;
        double d = 0.;
        double tol1, etemp;
        double CGOLD = 0.3819660;
        while (true) {
            double xm = 0.5*(a+b);
            double tol2 = 2.*(tol1=1.e-8*fabs(x)+1.e-10);
            if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
                return fx;
            }
            if (fabs(e) > tol1) {
                double r=(x-w)*(fx-fv);
                double q=(x-v)*(fx-fw);
                double p=(x-v)*q-(x-w)*r;
                q=2.*(q-r);
                if (q > 0.) p = -p;
                q=fabs(q);
                etemp = e;
                e = d;
                if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) 
                    d = CGOLD*(e=(x >= xm ? a-x : b-x));
                else {
                    d = p/q;
                    u = x+d;
                    if (u-a < tol2 || b-u < tol2)
                        d=copysign(tol1,xm-x);
                }
            } else {
                d = CGOLD*(e=(x >= xm ? a-x : b-x));
            }
            u = fabs(d) >= tol1 ? x+d : x+copysign(tol1,d);
            fu=distfunc(s,t,X(u),Y(u));
            if (fu <= fx) {
                if (u >= x) a=x; else b=x;
                v=w; w=x; x=u;
                fv=fw; fw=fx; fx=fu;
            } else {
                if (u < x) a=u; else b=u;
                if (fu <= fw || w == x) {
                    v=w;
                    w=u;
                    fv=fw;
                    fw=fu;
                } else if (fu <= fv || v == x || v == w) {
                    v=u;
                    fv=fu;
                }
            }
        }
    }

    double InitialFunc::XY(const double s, const double t) const
    {
        double dc = (stop-start)/xy_numsteps;
        double ang = 0;
        double val = DBL_MAX;
        double c;
        double x[2], y[2];
        double mx, my;
        int imin;
        double tval;

        x[0] = X(start);
        y[0] = Y(start);
        for (int i=0; i<xy_numsteps; ++i) {
            c = start+i*dc;
            x[1] = X(c+dc);
            y[1] = Y(c+dc);
            mx = (x[0]+x[1])/2;
            my = (y[0]+y[1])/2;
            tval = sqrt((s-mx)*(s-mx)+(t-my)*(t-my));
            if (tval < val) imin = i;
            val = min(val, tval);
            ang += asin(((y[1]-t)*(x[0]-s)-(y[0]-t)*(x[1]-s))/
                        sqrt((sqr(x[0]-s)+sqr(y[0]-t))
                             *(sqr(x[1]-s)+sqr(y[1]-t))));
            x[0] = x[1];
            y[0] = y[1];
        }
   
        val = Brent(start+imin*dc, dc, s, t);
      
        return val*(fabs(ang) > M_PI ? 1 : -1);
    }

}

#endif

      












