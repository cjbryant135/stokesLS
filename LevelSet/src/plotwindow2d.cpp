//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      PlotWindow2D Class
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
#include "plotwindow2d.h"
#ifndef M_PI
#define M_PI 3.1415926
#endif

namespace levelset {
        
//
//      ~PlotWindow2D
//
//      Denstructor for PlotWindow2D.  
//
//      Created 11/21/95
//
    PlotWindow2D::~PlotWindow2D(void)
    {
//      ~PlotWindow();
    }


//
//      PlotWindow2D::SetRange
//
//      Set the physical dimensions of the domain
//
//      Created 11/21/95
//
    void PlotWindow2D::SetRange(double xmin, double xmax, double ymin, double ymax)
    {
        pxmin = xmin;
        pxmax = xmax;
        pymin = ymin;
        pymax = ymax;
    }

//
//      PlotWindow2D::Line
//
//      Draw a line from (x0,y0) to (x1,y1)
//
//      Created 11/21/95
//
    void PlotWindow2D::Line(double x0, double y0, double x1, double y1)
    {
        PlotSegment(ispot(x0),jspot(y0),ispot(x1),jspot(y1));
    }

//
//      PlotWindow2D::Arrow
//
//      Draw an arrow from (x0,y0) pointing to (x1,y1)
//
//      Created 11/21/95
//
    void PlotWindow2D::Arrow(double x0, double y0, double x1, double y1)
    {
// arrow_frac is the fraction of the length the flanges should be.
        const double arrow_frac = 0.17;
// arrow_deg_angle is the angle of opening for the flanges in degrees.
        const double arrow_deg_angle = 30.;

        double angle = arrow_deg_angle*M_PI/180.;

        PlotSegment(ispot(x0),jspot(y0),ispot(x1),jspot(y1));
        double x2 = arrow_frac*(cos(angle)*(x0-x1)-sin(angle)*(y0-y1))+x1;
        double y2 = arrow_frac*(sin(angle)*(x0-x1)+cos(angle)*(y0-y1))+y1;
        PlotSegment(ispot(x1),jspot(y1),ispot(x2),jspot(y2));
        x2 = arrow_frac*(cos(angle)*(x0-x1)+sin(angle)*(y0-y1))+x1;
        y2 = arrow_frac*(-sin(angle)*(x0-x1)+cos(angle)*(y0-y1))+y1;
        PlotSegment(ispot(x1),jspot(y1),ispot(x2),jspot(y2));
    }

//
//      PlotWindow2D::Text
//
//      Print text at the physical coordinates
//
//      Created 11/21/95
//
    void PlotWindow2D::Text(double x0, double y0, const char *str)
    {
        PlotText(ispot(x0),jspot(y0),str);
    }

//
//      PlotWindow2D::Dot
//
//      Print a dot at the physical coordinates with the 
//      specified color (see PlotWindow.cp)
//
//      Created 11/21/95
//
    void PlotWindow2D::Dot(double x0, double y0, int color)
    {
        PlotDot(ispot(x0), jspot(y0), color);
    }

//
//      PlotWindow2D::Triangle
//
//      Plot a triangle with vertices at the given physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow2D::Triangle(double x0, double y0, double x1, double y1,
                                double x2, double y2, unsigned short color)
    {
        PlotTriangle(ispot(x0), jspot(y0), ispot(x1), jspot(y1), ispot(x2),
                     jspot(y2), color);
    }

//
//      PlotWindow2D::Polygon
//
//      Plot a polygon with vertices at the given physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow2D::Polygon(int num, double* x, double* y, unsigned short color)
    {
        int *ix = new int[num];
        int *iy = new int[num];
        for (int i=0; i<num; ++i) {
            ix[i] = ispot(x[i]);
            iy[i] = jspot(y[i]);
        }
        PlotPolygon(num, ix, iy, color);
        delete[] ix;
        delete[] iy;
    }

//
//      PlotWindow2D::Circle
//
//      Plot a circle centered at (x0,y0) with radius r in the physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow2D::Circle(double x0, double y0, double radius)
    {
        PlotOval(ispot(x0), jspot(y0), ispot(x0+radius)-ispot(x0), 
                 jspot(y0+radius)-jspot(y0));
    }

//
//      PlotWindow2D::FilledCircle
//
//      Plot a circle centered at (x0,y0) with radius r in the physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow2D::FilledCircle(double x0, double y0, double radius)
    {
        PlotFilledOval(ispot(x0), jspot(y0), ispot(x0+radius)-ispot(x0), 
                       jspot(y0+radius)-jspot(y0));
    }

//
//      PlotWindow2D::Clear
//
//      Clear a rectangle at (x0,y0) with width w, and height h
//
//      Created 8/10/04
//
    void PlotWindow2D::Clear(double x0, double y0, double w, double h)
    {
        PlotWindow::Clear(ispot(x0), jspot(y0), ispot(x0+w)-ispot(x0), jspot(y0+h)-jspot(y0));
    }

}
