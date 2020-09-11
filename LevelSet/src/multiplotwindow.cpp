//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      MultiPlotWindow2D Class
//
//      This class is designed to handle the background operations of
//      running a graphics window.  To take advantage of the primitives,
//      a subclass should be written to handle the interface to the
//      application.
//
//      Created 11/21/95
//
// Modified 4/11/97 - Add X Windows functionality
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <math.h>
#include <limits.h>
#include <float.h>
#include "multiplotwindow.h"
#include "utility.h"
#ifndef M_PI
#define M_PI 3.1415926
#endif
typedef bool Boolean;

namespace levelset {
        
    MultiPlotWindow::MultiPlotWindow(const int r, const int c, Boolean yes_window, Boolean pae,
                                     const int w, const int h)
        : mpxmin(NULL), mpymin(NULL), mpxmax(NULL), mpymax(NULL), rows(r), cols(c),
          cur_r(0), cur_c(0), PlotWindow3D(yes_window, pae, w, h)
    {
        mpxmin = new double*[rows];
        mpxmax = new double*[rows];
        mpymin = new double*[rows];
        mpymax = new double*[rows];
        for (int i=0; i<rows; ++i) {
            mpxmin[i] = new double[cols];
            mpxmax[i] = new double[cols];
            mpymin[i] = new double[cols];
            mpymax[i] = new double[cols];
        }
        for (int i=0; i<3; ++i) {
            for (int j=0; j<2; ++j) {
                limits[i][j] = new double*[rows];
                for (int k=0; k<rows; ++k)
                    limits[i][j][k] = new double[cols];
            }
            view[i] = new double*[rows];
            for (int k=0; k<rows; ++k)
                view[i][k] = new double[cols];
            b_[i] = new double*[rows];
            for (int k=0; k<rows; ++k)
                b_[i][k] = new double[cols];
        }
        for (int i=0; i<2; ++i) {
            a_[i] = new double*[rows];
            for (int k=0; k<rows; ++k)
                a_[i][k] = new double[cols];
        }
        vdv = new double*[rows];
        for (int i=0; i<rows; ++i)
            vdv[i] = new double[cols];
        PlotWindow2D::SetRange(0.,double(cols),0.,double(rows));
        for (int i=1; i<rows; ++i) 
            PlotWindow2D::Line(0,double(i),cols,double(i));
        for (int i=1; i<cols; ++i)
            PlotWindow2D::Line(double(i),0,double(i),rows);
    }

    MultiPlotWindow::~MultiPlotWindow(void)
    {
        if (mpxmin != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpxmin[i];
            delete[] mpxmin;
        }
        if (mpxmax != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpxmax[i];
            delete[] mpxmax;
        }
        if (mpymin != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpymin[i];
            delete[] mpymin;
        }
        if (mpymax != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpymax[i];
            delete[] mpymax;
        }
        if (limits[0][0] != NULL) {
            for (int i=0; i<3; ++i)
                for (int j=0; j<2; ++j) {
                    for (int k=0; k<rows; ++k)
                        delete[] limits[i][j][k];
                    delete[] limits[i][j];
                }
        }
        if (view[0] != NULL) {
            for (int i=0; i<3; ++i) {
                for (int j=0; j<rows; ++j)
                    delete[] view[i][j];
                delete[] view[i];
            }
        }
        if (vdv != NULL) {
            for (int i=0; i<rows; ++i)
                delete[] vdv[i];
            delete[] vdv;
        }
        if (a_[0] != NULL) {
            for (int i=0; i<2; ++i) {
                for (int j=0; j<rows; ++j)
                    delete[] a_[i][j];
                delete[] a_[i];
            }
        }
        if (b_[0] != NULL) {
            for (int i=0; i<3; ++i) {
                for (int j=0; j<rows; ++j)
                    delete[] b_[i][j];
                delete[] b_[i];
            }
        }
    }

    void MultiPlotWindow::LocalToGlobal(const int r, const int c, const double x,
                                        const double y, double& a, double& b) const
    {
        a = c+(x-mpxmin[r][c])/(mpxmax[r][c]-mpxmin[r][c]);
        b = r+(y-mpymin[r][c])/(mpymax[r][c]-mpymin[r][c]);
    }

    void MultiPlotWindow::SetRange(const int r, const int c, double xmin, double xmax, double ymin, double ymax)
    {
        mpxmin[r][c] = xmin;
        mpxmax[r][c] = xmax;
        mpymin[r][c] = ymin;
        mpymax[r][c] = ymax;
    }

    void MultiPlotWindow::Line(const int r, const int c, double x0, double y0, double x1, double y1)
    {
        double a0, b0, a1, b1;
        LocalToGlobal(r,c,x0,y0,a0,b0);
        LocalToGlobal(r,c,x1,y1,a1,b1);
        PlotWindow2D::Line(a0,b0,a1,b1);
    }

    void MultiPlotWindow::Arrow(const int r, const int c, double x0, double y0, double x1, double y1)
    {
        double a0, b0, a1, b1;
        LocalToGlobal(r,c,x0,y0,a0,b0);
        LocalToGlobal(r,c,x1,y1,a1,b1);
        PlotWindow2D::Arrow(a0,b0,a1,b1);
    }

    void MultiPlotWindow::Text(const int r, const int c, double x0, double y0, const char *str)
    {
        double a, b;
        LocalToGlobal(r,c,x0,y0,a,b);
        PlotWindow2D::Text(a,b,str);
    }

    void MultiPlotWindow::Text(const int r, const int c, PlotWindow::WinPos wp, const char* s)
    {
        switch(wp) {
        case LL: PlotWindow2D::Text(c+0.02, r+0.02, s); break;
        case LR: PlotWindow2D::Text(c+0.9, r+0.02, s); break;
        case UL: PlotWindow2D::Text(c+0.02, r+0.95, s); break;
        case UR: PlotWindow2D::Text(c+0.9, r+0.95, s); break;
        }
    }

    void MultiPlotWindow::Dot(const int r, const int c, double x0, double y0, int color)
    {
        double a, b;
        LocalToGlobal(r,c,x0,y0,a,b);
        PlotWindow2D::Dot(a, b, color);
    }

    void MultiPlotWindow::Triangle(const int r, const int c, double x0, double y0, double x1, double y1,
                                   double x2, double y2, unsigned short color)
    {
        double a[3], b[3];
        LocalToGlobal(r,c,x0,y0,a[0],b[0]);
        LocalToGlobal(r,c,x1,y1,a[1],b[1]);
        LocalToGlobal(r,c,x2,y2,a[2],b[2]);
        PlotWindow2D::Triangle(a[0],b[0],a[1],b[1],a[2],b[2],color);
    }

    void MultiPlotWindow::Polygon(const int r, const int c, int num, double* x, double* y, unsigned short color)
    {
        double *a = new double[num];
        double *b = new double[num];
        for (int i=0; i<num; ++i) 
            LocalToGlobal(r,c,x[i],y[i],a[i],b[i]);
        PlotWindow2D::Polygon(num, a, b, color);
        delete[] a;
        delete[] b;
    }

    void MultiPlotWindow::Circle(const int r, const int c, double x0, double y0, double radius)
    {
        double a, b, d;
        LocalToGlobal(r,c,x0,y0,a,b);
        LocalToGlobal(r,c,x0+radius,y0,d,b);
        d = d-a;
        PlotWindow2D::Circle(a,b,d);
    }

    void MultiPlotWindow::Clear(const int r, const int c, double x0, double y0, double w, double h)
    {
        double a, b, ww, hh;
        LocalToGlobal(r,c,x0,y0,a,b);
        LocalToGlobal(r,c,x0+w,y0+h,ww,hh);
        ww = ww - a;
        hh = hh - b;
        PlotWindow2D::Clear(a,b,ww,hh);
    }

    void MultiPlotWindow::Clear(const int r, const int c)
    {
        PlotWindow2D::Clear(c,r,c+1,r+1);
        PlotWindow2D::Line(c,r,c+1,r);
        PlotWindow2D::Line(c+1,r,c+1,r+1);
        PlotWindow2D::Line(c+1,r+1,c,r+1);
        PlotWindow2D::Line(c,r+1,c,r);
    }

    void MultiPlotWindow::SetRange(const int r, const int c, double xmin, double xmax,
                                   double ymin, double ymax, double zmin, double zmax)
    {
        limits[0][0][r][c] = xmin;
        limits[0][1][r][c] = xmax;
        limits[1][0][r][c] = ymin;
        limits[1][1][r][c] = ymax;
        limits[2][0][r][c] = zmin;
        limits[2][1][r][c] = zmax;
        Recalc(r,c);
    }

    void MultiPlotWindow::GetRange(const int r, const int c, double& xmin, double& xmax,
                                   double& ymin, double& ymax, double& zmin, double& zmax) const
    {
        xmin = limits[0][0][r][c];
        xmax = limits[0][1][r][c];
        ymin = limits[1][0][r][c];
        ymax = limits[1][1][r][c];
        zmin = limits[2][0][r][c];
        zmax = limits[2][1][r][c];
    }

    void MultiPlotWindow::GetActualRange(const int r, const int c, double& xmin, double& xmax, double& ymin,
                                         double& ymax, double& zmin, double& zmax) const
    {
        xmin = -DBL_MAX;
        xmax = DBL_MAX;
        ymin = -DBL_MAX;
        ymax = DBL_MAX;
        zmin = -DBL_MAX;
        zmax = DBL_MAX;
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j) {
                double vmax = (mpxmax[r][c]*(vdv[r][c]-limits[1][i][r][c]*view[1][r][c]-limits[2][j][r][c]*view[2][r][c])
                               -a_[1][r][c]*limits[1][i][r][c])/(mpxmax[r][c]*view[0][r][c]+a_[0][r][c]);
                double vmin = (mpxmin[r][c]*(vdv[r][c]-limits[1][i][r][c]*view[1][r][c]-limits[2][j][r][c]*view[2][r][c])
                               -a_[1][r][c]*limits[1][i][r][c])/(mpxmin[r][c]*view[0][r][c]+a_[0][r][c]);
                double wmax = (mpymax[r][c]*(vdv[r][c]-limits[1][i][r][c]*view[1][r][c]-limits[2][j][r][c]*view[2][r][c])
                               -b_[1][r][c]*limits[1][i][r][c]-b_[2][r][c]*limits[2][j][r][c])/(mpymax[r][c]*view[0][r][c]+b_[0][r][c]);
                double wmin = (mpymin[r][c]*(vdv[r][c]-limits[1][i][r][c]*view[1][r][c]-limits[2][j][r][c]*view[2][r][c])
                               -b_[1][r][c]*limits[1][i][r][c]-b_[2][r][c]*limits[2][j][r][c])/(mpymin[r][c]*view[0][r][c]+b_[0][r][c]);
                xmax = min(xmax,max(vmax,vmin),max(wmax,wmin));
                xmin = max(xmin,min(vmax,vmin),min(wmax,wmin));
                vmax = (mpxmax[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[2][j][r][c]*view[2][r][c])
                        -a_[0][r][c]*limits[0][i][r][c])/(mpxmax[r][c]*view[1][r][c]+a_[1][r][c]);
                vmin = (mpxmin[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[2][j][r][c]*view[2][r][c])
                        -a_[0][r][c]*limits[0][i][r][c])/(mpxmin[r][c]*view[1][r][c]+a_[1][r][c]);
                wmax = (mpymax[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[2][j][r][c]*view[2][r][c])
                        -b_[0][r][c]*limits[0][i][r][c]-b_[2][r][c]*limits[2][j][r][c])/(mpymax[r][c]*view[1][r][c]+b_[1][r][c]);
                wmin = (mpymin[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[2][j][r][c]*view[2][r][c])
                        -b_[0][r][c]*limits[0][i][r][c]-b_[2][r][c]*limits[2][j][r][c])/(mpymin[r][c]*view[1][r][c]+b_[1][r][c]);
                ymax = min(ymax,max(vmax,vmin),max(wmax,wmin));
                ymin = max(ymin,min(vmax,vmin),min(wmax,wmin));
#if 0
                vmax = (mpxmax[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[1][j][r][c]*view[1][r][c])
                        -a_[0][r][c]*limits[0][i][r][c]-a_[1][r][c]*limits[1][j][r][c])/(mpxmax[r][c]*view[2][r][c]);
                vmin = (mpxmin[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[1][j][r][c]*view[1][r][c])
                        -a_[0][r][c]*limits[0][i][r][c]-a_[1][r][c]*limits[1][j][r][c])/(mpxmin[r][c]*view[2][r][c]);
#endif
                wmax = (mpymax[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[1][j][r][c]*view[1][r][c])
                        -b_[0][r][c]*limits[0][i][r][c]-b_[1][r][c]*limits[1][j][r][c])/(mpymax[r][c]*view[2][r][c]+b_[2][r][c]);
                wmin = (mpymin[r][c]*(vdv[r][c]-limits[0][i][r][c]*view[0][r][c]-limits[1][j][r][c]*view[1][r][c])
                        -b_[0][r][c]*limits[0][i][r][c]-b_[1][r][c]*limits[1][j][r][c])/(mpymin[r][c]*view[2][r][c]+b_[2][r][c]);
                zmax = min(zmax,max(wmax,wmin));
                zmin = max(zmin,min(wmax,wmin));
            }
    }                   

    void MultiPlotWindow::SetView(const int r, const int c, double x, double y, double z)
    {
        view[0][r][c] = x;
        view[1][r][c] = y;
        view[2][r][c] = z;
        Recalc(r,c);
    }

    void MultiPlotWindow::GetView(const int r, const int c, double& x, double& y, double& z) const
    {
        x = view[0][r][c];
        y = view[1][r][c];
        z = view[2][r][c];
    }

    void MultiPlotWindow::Axes(const int r, const int c)
    {
        Line(r,c,limits[0][0][r][c],limits[1][0][r][c],limits[2][0][r][c],
             limits[0][1][r][c],limits[1][0][r][c],limits[2][0][r][c]);
        Text(r,c,limits[0][1][r][c],limits[1][0][r][c],limits[2][0][r][c],"X");
        Line(r,c,limits[0][0][r][c],limits[1][0][r][c],limits[2][0][r][c],
             limits[0][0][r][c],limits[1][1][r][c],limits[2][0][r][c]);
        Text(r,c,limits[0][0][r][c],limits[1][1][r][c],limits[2][0][r][c],"Y");
        Line(r,c,limits[0][0][r][c],limits[1][0][r][c],limits[2][0][r][c],
             limits[0][0][r][c],limits[1][0][r][c],limits[2][1][r][c]);
        Text(r,c,limits[0][0][r][c],limits[1][0][r][c],limits[2][1][r][c],"Z");
    }

    void MultiPlotWindow::Line(const int r, const int c, double x0, double y0, double z0,
                               double x1, double y1, double z1)
    {
        Line(r,c,projx(r,c,x0,y0,z0),projy(r,c,x0,y0,z0),
             projx(r,c,x1,y1,z1),projy(r,c,x1,y1,z1));
    }

    void MultiPlotWindow::Arrow(const int r, const int c, double x0, double y0, double z0,
                                double x1, double y1, double z1)
    {
        Arrow(r,c,projx(r,c,x0,y0,z0),projy(r,c,x0,y0,z0),
              projx(r,c,x1,y1,z1),projy(r,c,x1,y1,z1));
    }

    void MultiPlotWindow::Text(const int r, const int c, double x0, double y0, double z0,
                               char *str)
    {
        Text(r,c,projx(r,c,x0,y0,z0),projy(r,c,x0,y0,z0),str);
    }

    void MultiPlotWindow::Dot(const int r, const int c, double x0, double y0, double z0, 
                              int color)
    {
        Dot(r,c,projx(r,c,x0,y0,z0),projy(r,c,x0,y0,z0),color);
    }

    void MultiPlotWindow::Triangle(const int r, const int c, double x0, double y0, 
                                   double z0, double x1, double y1, double z1, double x2, double y2, 
                                   double z2)
    {
        double vec[3];
        double norm;
        unsigned short shade;
        
        vec[0] = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
        vec[1] = (x2-x0)*(z1-z0)-(x1-x0)*(z2-z0);
        vec[2] = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
        norm = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
        shade = (unsigned short)(USHRT_MAX*fabs((vec[0]*view[0][r][c]+vec[1]*view[1][r][c]
                                                 +vec[2]*view[2][r][c])/norm/sqrt(vdv[r][c])));
        Triangle(r,c,projx(r,c,x0,y0,z0), projy(r,c,x0,y0,z0),
                 projx(r,c,x1,y1,z1), projy(r,c,x1,y1,z1),
                 projx(r,c,x2,y2,z2), projy(r,c,x2,y2,z2), shade);
    }

    void MultiPlotWindow::Polygon(const int r, const int c, int num, 
                                  double* x, double* y, double* z, unsigned short color)
    {
        double *tx = new double[num];
        double *ty = new double[num];
        for (int i=0; i<num; ++i) {
            tx[i] = projx(r,c,x[i],y[i],z[i]);
            ty[i] = projy(r,c,x[i],y[i],z[i]);
        }
        Polygon(r, c, num, tx, ty, color);
        delete[] tx;
        delete[] ty;
    }

    void MultiPlotWindow::Circle(const int r, const int c, double x0, double y0,
                                 double z0, double radius)
    {
        Circle(r,c,projx(r,c,x0,y0,z0),projy(r,c,x0,y0,z0), radius);
    }

    void MultiPlotWindow::Recalc(const int r, const int c)
    {   
        double svdv, denom;
   
        vdv[r][c] = view[0][r][c]*view[0][r][c]+view[1][r][c]*view[1][r][c]+view[2][r][c]*view[2][r][c];
        svdv = sqrt(vdv[r][c]);
        denom = sqrt(view[0][r][c]*view[0][r][c]+view[1][r][c]*view[1][r][c]);
        a_[0][r][c] = -vdv[r][c]*view[1][r][c]/denom;
        a_[1][r][c] = vdv[r][c]*view[0][r][c]/denom;
        b_[0][r][c] = -svdv*view[0][r][c]*view[2][r][c]/denom;
        b_[1][r][c] = -svdv*view[1][r][c]*view[2][r][c]/denom;
        b_[2][r][c] = svdv*denom;

        double ptemp[8][2];
        int n = 0;
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j)
                for (int k=0; k<2; ++k) {
                    ptemp[n][0] = projx(r,c,limits[0][i][r][c], limits[1][j][r][c], limits[2][k][r][c]);
                    ptemp[n++][1] = projy(r,c,limits[0][i][r][c], limits[1][j][r][c], limits[2][k][r][c]);
                }
        mpxmax[r][c] = max(max(ptemp[0][0],ptemp[1][0],ptemp[2][0],ptemp[3][0]),
                           max(ptemp[4][0],ptemp[5][0],ptemp[6][0],ptemp[7][0]));
        mpxmin[r][c] = min(min(ptemp[0][0],ptemp[1][0],ptemp[2][0],ptemp[3][0]),
                           min(ptemp[4][0],ptemp[5][0],ptemp[6][0],ptemp[7][0]));
        mpymax[r][c] = max(max(ptemp[0][1],ptemp[1][1],ptemp[2][1],ptemp[3][1]),
                           max(ptemp[4][1],ptemp[5][1],ptemp[6][1],ptemp[7][1]));
        mpymin[r][c] = min(min(ptemp[0][1],ptemp[1][1],ptemp[2][1],ptemp[3][1]),
                           min(ptemp[4][1],ptemp[5][1],ptemp[6][1],ptemp[7][1]));

        double xdiff = mpxmax[r][c] - mpxmin[r][c];
        double ydiff = mpymax[r][c] - mpymin[r][c];
        double temp = max(xdiff, ydiff);
        mpxmax[r][c] = mpxmax[r][c] + (temp-xdiff)/2.;
        mpxmin[r][c] = mpxmin[r][c] - (temp-xdiff)/2.;
        mpymax[r][c] = mpymax[r][c] + (temp-ydiff)/2.;
        mpymin[r][c] = mpymin[r][c] - (temp-ydiff)/2.;
    }

}

