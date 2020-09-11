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
#include "multiplotwindow2d.h"
#ifndef M_PI
#define M_PI 3.1415926
#endif
typedef bool Boolean;

namespace levelset {
        
    MultiPlotWindow2D::MultiPlotWindow2D(const int r, const int c, Boolean yes_window, Boolean pae,
                                         const int w, const int h)
        : mpxmin(NULL), mpymin(NULL), mpxmax(NULL), mpymax(NULL), rows(r), cols(c),
          cur_r(0), cur_c(0), PlotWindow2D(yes_window, pae, w, h)
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
        PlotWindow2D::SetRange(0.,double(cols),0.,double(rows));
        for (int i=1; i<rows; ++i) 
            PlotWindow2D::Line(0,double(i),cols,double(i));
        for (int i=1; i<cols; ++i)
            PlotWindow2D::Line(double(i),0,double(i),rows);
    }

    MultiPlotWindow2D::~MultiPlotWindow2D(void)
    {
        if (mpxmin != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpxmin[0];
            delete[] mpxmin;
        }
        if (mpxmax != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpxmax[0];
            delete[] mpxmax;
        }
        if (mpymin != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpymin[0];
            delete[] mpymin;
        }
        if (mpymax != NULL) {
            for (int i=0; i<rows; ++i) 
                delete[] mpymax[0];
            delete[] mpymax;
        }
    }

    void MultiPlotWindow2D::LocalToGlobal(const int r, const int c, const double x,
                                          const double y, double& a, double& b) const
    {
        a = c+(x-mpxmin[r][c])/(mpxmax[r][c]-mpxmin[r][c]);
        b = r+(y-mpymin[r][c])/(mpymax[r][c]-mpymin[r][c]);
    }

    void MultiPlotWindow2D::SetRange(const int r, const int c, double xmin, double xmax, double ymin, double ymax)
    {
        mpxmin[r][c] = xmin;
        mpxmax[r][c] = xmax;
        mpymin[r][c] = ymin;
        mpymax[r][c] = ymax;
    }

    void MultiPlotWindow2D::Line(const int r, const int c, double x0, double y0, double x1, double y1)
    {
        double a0, b0, a1, b1;
        LocalToGlobal(r,c,x0,y0,a0,b0);
        LocalToGlobal(r,c,x1,y1,a1,b1);
        PlotWindow2D::Line(a0,b0,a1,b1);
    }

    void MultiPlotWindow2D::Arrow(const int r, const int c, double x0, double y0, double x1, double y1)
    {
        double a0, b0, a1, b1;
        LocalToGlobal(r,c,x0,y0,a0,b0);
        LocalToGlobal(r,c,x1,y1,a1,b1);
        PlotWindow2D::Arrow(a0,b0,a1,b1);
    }

    void MultiPlotWindow2D::Text(const int r, const int c, double x0, double y0, const char *str)
    {
        double a, b;
        LocalToGlobal(r,c,x0,y0,a,b);
        PlotWindow2D::Text(a,b,str);
    }

    void MultiPlotWindow2D::Text(const int r, const int c, PlotWindow::WinPos wp, const char* s)
    {
        switch(wp) {
        case LL: PlotWindow2D::Text(c+0.02, r+0.02, s); break;
        case LR: PlotWindow2D::Text(c+0.9, r+0.02, s); break;
        case UL: PlotWindow2D::Text(c+0.02, r+0.95, s); break;
        case UR: PlotWindow2D::Text(c+0.9, r+0.95, s); break;
        }
    }

    void MultiPlotWindow2D::Dot(const int r, const int c, double x0, double y0, int color)
    {
        double a, b;
        LocalToGlobal(r,c,x0,y0,a,b);
        PlotWindow2D::Dot(a, b, color);
    }

    void MultiPlotWindow2D::Triangle(const int r, const int c, double x0, double y0, double x1, double y1,
                                     double x2, double y2, unsigned short color)
    {
        double a[3], b[3];
        LocalToGlobal(r,c,x0,y0,a[0],b[0]);
        LocalToGlobal(r,c,x1,y1,a[1],b[1]);
        LocalToGlobal(r,c,x2,y2,a[2],b[2]);
        PlotWindow2D::Triangle(a[0],b[0],a[1],b[1],a[2],b[2],color);
    }

    void MultiPlotWindow2D::Polygon(const int r, const int c, int num, double* x, double* y, unsigned short color)
    {
        double *a = new double[num];
        double *b = new double[num];
        for (int i=0; i<num; ++i) 
            LocalToGlobal(r,c,x[i],y[i],a[i],b[i]);
        PlotWindow2D::Polygon(num, a, b, color);
        delete[] a;
        delete[] b;
    }

    void MultiPlotWindow2D::Circle(const int r, const int c, double x0, double y0, double radius)
    {
        double a, b, d;
        LocalToGlobal(r,c,x0,y0,a,b);
        LocalToGlobal(r,c,x0+radius,y0,d,b);
        d = d-a;
        PlotWindow2D::Circle(a,b,d);
    }

    void MultiPlotWindow2D::Clear(const int r, const int c, double x0, double y0, double w, double h)
    {
        double a, b, ww, hh;
        LocalToGlobal(r,c,x0,y0,a,b);
        LocalToGlobal(r,c,x0+w,y0+h,ww,hh);
        ww = ww - a;
        hh = hh - b;
        PlotWindow2D::Clear(a,b,ww,hh);
    }

    void MultiPlotWindow2D::Clear(const int r, const int c)
    {
        PlotWindow2D::Clear(c,c+1,r,r+1);
        PlotWindow2D::Line(c,r,c+1,r);
        PlotWindow2D::Line(c+1,r,c+1,r+1);
        PlotWindow2D::Line(c+1,r+1,c,r+1);
        PlotWindow2D::Line(c,r+1,c,r);
    }

}
