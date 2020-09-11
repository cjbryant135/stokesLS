#include "mound.h"
#include "utility.h"
#include <float.h>

namespace levelset {
        
    void Mound::SetParams(InputParams *param)
    {
        char fname[256];
        param->GetCharParam("Curve Data File", fname, "");
        ReadData(fname);
    }

    void Mound::CleanUp(void)
    {
        int i;
        
        if (x != NULL) {
            for (i=0; i<numc; ++i) delete[] x[i];
            delete[] x;
        }
        if (y != NULL) {
            for (i=0; i<numc; ++i) delete[] y[i];
            delete[] y;
        }
        if (num != NULL) 
            delete[] num;
    }

    void Mound::ReadData(const char* fname)
    {
        std::ifstream data(fname);
        data >> numc;
        if (num != NULL) delete[] num;
        num = new int[numc];
        if (x != NULL) delete[] x;
        x = new double*[numc];
        if (y != NULL) delete[] y;
        y = new double*[numc];
        for (int i=0; i<numc; ++i) {
            data >> num[i];
            x[i] = new double[num[i]+1];
            y[i] = new double[num[i]+1];
            for (int j=0; j<num[i]; ++j) 
                data >> x[i][j] >> y[i][j];
            x[i][num[i]] = x[i][0];
            y[i][num[i]] = y[i][0];
        }
    }

    double dist2seg(const double sx1, const double sy1, const double sx2, const double sy2, 
                    const double x, const double y)
    {
        double dx[2], dy[2];
        dx[0] = sx2-sx1;
        dx[1] = x-sx1;
        dy[0] = sy2-sy1;
        dy[1] = y-sy1;
        double t = (dx[0]*dx[1]+dy[0]*dy[1])/(sqr(dx[0])+sqr(dy[0]));
        if (t < 0.) 
            return sqrt(sqr(dx[1])+sqr(dy[1]));
        else if (t > 1.) 
            return sqrt(sqr(x-sx2)+sqr(y-sy2));
        else
            return sqrt((sqr(dx[0])*sqr(dy[1])+sqr(dy[0])*sqr(dx[1])
                         -2*dx[0]*dx[1]*dy[0]*dy[1])/(sqr(dx[0])+sqr(dy[0])));
    }

    double Mound::XYZ(const double s, const double t, const double u) const
    {
        int i, j;
        double ans = DBL_MAX;
        double d;
        double tans = -DBL_MAX;
        int sgn = -1;
        int tsgn = -1;
        for (i=0; i<numc; ++i) {
            ans = DBL_MAX;
            sgn = -1;
            for (j=0; j<num[i]; ++j) {
                d = dist2seg(x[i][j],y[i][j],x[i][j+1],y[i][j+1],s,t);
                ans = ans < d ? ans : d;
                if ((s == x[i][j] || (s-x[i][j])*(s-x[i][j+1]) < 0)
                    && t > y[i][j]+(y[i][j+1]-y[i][j])/(x[i][j+1]-x[i][j])*(s-x[i][j])) 
                    sgn = -sgn;
            }
            tsgn = max(tsgn, sgn);
            tans = max(tans, sgn*ans);
        }
        
        return tans-u;
    }

}

