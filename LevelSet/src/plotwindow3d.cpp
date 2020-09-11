//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      PlotWindow3D Class
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
#include "plotwindow3d.h"
#include "utility.h"
#ifndef M_PI
#define M_PI 3.1415926
#endif

namespace levelset {
        
//
//      ~PlotWindow3D
//
//      Denstructor for PlotWindow3D.  
//
//      Created 11/21/95
//
    PlotWindow3D::~PlotWindow3D(void)
    {
//      ~PlotWindow();
    }


//
//      PlotWindow3D::SetRange
//
//      Set the physical dimensions of the domain
//
//      Created 11/21/95
//
    void PlotWindow3D::SetRange(double xmin, double xmax, double ymin, double ymax,
                                double zmin, double zmax)
    {
        limits[0][0] = xmin;
        limits[0][1] = xmax;
        limits[1][0] = ymin;
        limits[1][1] = ymax;
        limits[2][0] = zmin;
        limits[2][1] = zmax;
        Recalc();
    }

    void PlotWindow3D::SetView(double x, double y, double z)
    {
        view[0] = x;
        view[1] = y;
        view[2] = z;
        Recalc();
    }

    void PlotWindow3D::GetRange(double& xmin, double& xmax, double& ymin,
                                double& ymax, double& zmin, double& zmax) const
    {
        xmin = limits[0][0];
        xmax = limits[0][1];
        ymin = limits[1][0];
        ymax = limits[1][1];
        zmin = limits[2][0];
        zmax = limits[2][1];
    }

    void PlotWindow3D::GetActualRange(double& xmin, double& xmax, double& ymin,
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
                double vmax = (pxmax*(vdv-limits[1][i]*view[1]-limits[2][j]*view[2])
                               -a_[1]*limits[1][i])/(pxmax*view[0]+a_[0]);
                double vmin = (pxmin*(vdv-limits[1][i]*view[1]-limits[2][j]*view[2])
                               -a_[1]*limits[1][i])/(pxmin*view[0]+a_[0]);
                double wmax = (pymax*(vdv-limits[1][i]*view[1]-limits[2][j]*view[2])
                               -b_[1]*limits[1][i]-b_[2]*limits[2][j])/(pymax*view[0]+b_[0]);
                double wmin = (pymin*(vdv-limits[1][i]*view[1]-limits[2][j]*view[2])
                               -b_[1]*limits[1][i]-b_[2]*limits[2][j])/(pymin*view[0]+b_[0]);
                xmax = min(xmax,max(vmax,vmin),max(wmax,wmin));
                xmin = max(xmin,min(vmax,vmin),min(wmax,wmin));
                vmax = (pxmax*(vdv-limits[0][i]*view[0]-limits[2][j]*view[2])
                        -a_[0]*limits[0][i])/(pxmax*view[1]+a_[1]);
                vmin = (pxmin*(vdv-limits[0][i]*view[0]-limits[2][j]*view[2])
                        -a_[0]*limits[0][i])/(pxmin*view[1]+a_[1]);
                wmax = (pymax*(vdv-limits[0][i]*view[0]-limits[2][j]*view[2])
                        -b_[0]*limits[0][i]-b_[2]*limits[2][j])/(pymax*view[1]+b_[1]);
                wmin = (pymin*(vdv-limits[0][i]*view[0]-limits[2][j]*view[2])
                        -b_[0]*limits[0][i]-b_[2]*limits[2][j])/(pymin*view[1]+b_[1]);
                ymax = min(ymax,max(vmax,vmin),max(wmax,wmin));
                ymin = max(ymin,min(vmax,vmin),min(wmax,wmin));
//                      vmax = (pxmax*(vdv-limits[0][i]*view[0]-limits[1][j]*view[1])
//                                              -a_[0]*limits[0][i]-a_[1]*limits[1][j])/(pxmax*view[2]);
//                      vmin = (pxmin*(vdv-limits[0][i]*view[0]-limits[1][j]*view[1])
//                                              -a_[0]*limits[0][i]-a_[1]*limits[1][j])/(pxmin*view[2]);
                wmax = (pymax*(vdv-limits[0][i]*view[0]-limits[1][j]*view[1])
                        -b_[0]*limits[0][i]-b_[1]*limits[1][j])/(pymax*view[2]+b_[2]);
                wmin = (pymin*(vdv-limits[0][i]*view[0]-limits[1][j]*view[1])
                        -b_[0]*limits[0][i]-b_[1]*limits[1][j])/(pymin*view[2]+b_[2]);
                zmax = min(zmax,max(wmax,wmin));
                zmin = max(zmin,min(wmax,wmin));
            }
    }                   

    void PlotWindow3D::GetView(double& x, double& y, double& z) const
    {
        x = view[0];
        y = view[1];
        z = view[2];
    }

    void PlotWindow3D::Axes(void) 
    {
        Line(limits[0][0],limits[1][0],limits[2][0],
             limits[0][1],limits[1][0],limits[2][0]);
        Text(limits[0][1],limits[1][0],limits[2][0],"X");
        Line(limits[0][0],limits[1][0],limits[2][0],
             limits[0][0],limits[1][1],limits[2][0]);
        Text(limits[0][0],limits[1][1],limits[2][0],"Y");
        Line(limits[0][0],limits[1][0],limits[2][0],
             limits[0][0],limits[1][0],limits[2][1]);
        Text(limits[0][0],limits[1][0],limits[2][1],"Z");
    }

    void PlotWindow3D::Recalc(void)
    {   
        double svdv, denom;
   
        vdv = view[0]*view[0]+view[1]*view[1]+view[2]*view[2];
        svdv = sqrt(vdv);
        denom = sqrt(view[0]*view[0]+view[1]*view[1]);
        a_[0] = -vdv*view[1]/denom;
        a_[1] = vdv*view[0]/denom;
        b_[0] = -svdv*view[0]*view[2]/denom;
        b_[1] = -svdv*view[1]*view[2]/denom;
        b_[2] = svdv*denom;

        double ptemp[8][2];
        int n = 0;
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j)
                for (int k=0; k<2; ++k) {
                    ptemp[n][0] = projx(limits[0][i], limits[1][j], limits[2][k]);
                    ptemp[n++][1] = projy(limits[0][i], limits[1][j], limits[2][k]);
                }
        pxmax = max(max(ptemp[0][0],ptemp[1][0],ptemp[2][0],ptemp[3][0]),
                    max(ptemp[4][0],ptemp[5][0],ptemp[6][0],ptemp[7][0]));
        pxmin = min(min(ptemp[0][0],ptemp[1][0],ptemp[2][0],ptemp[3][0]),
                    min(ptemp[4][0],ptemp[5][0],ptemp[6][0],ptemp[7][0]));
        pymax = max(max(ptemp[0][1],ptemp[1][1],ptemp[2][1],ptemp[3][1]),
                    max(ptemp[4][1],ptemp[5][1],ptemp[6][1],ptemp[7][1]));
        pymin = min(min(ptemp[0][1],ptemp[1][1],ptemp[2][1],ptemp[3][1]),
                    min(ptemp[4][1],ptemp[5][1],ptemp[6][1],ptemp[7][1]));

        double xdiff = pxmax - pxmin;
        double ydiff = pymax - pymin;
        double temp = max(xdiff, ydiff);
        pxmax = pxmax + (temp-xdiff)/2.;
        pxmin = pxmin - (temp-xdiff)/2.;
        pymax = pymax + (temp-ydiff)/2.;
        pymin = pymin - (temp-ydiff)/2.;
    }

//
//      PlotWindow3D::Line
//
//      Draw a line from (x0,y0) to (x1,y1)
//
//      Created 11/21/95
//
    void PlotWindow3D::Line(double x0, double y0, double z0,
                            double x1, double y1, double z1)
    {
        PlotWindow2D::Line(projx(x0,y0,z0),projy(x0,y0,z0),
                           projx(x1,y1,z1),projy(x1,y1,z1));
    }

//
//      PlotWindow3D::Arrow
//
//      Draw an arrow from (x0,y0) pointing to (x1,y1)
//
//      Created 11/21/95
//
    void PlotWindow3D::Arrow(double x0, double y0, double z0,
                             double x1, double y1, double z1)
    {
        PlotWindow2D::Arrow(projx(x0,y0,z0),projy(x0,y0,z0),
                            projx(x1,y1,z1),projy(x1,y1,z1));
    }

//
//      PlotWindow3D::Text
//
//      Print text at the physical coordinates
//
//      Created 11/21/95
//
    void PlotWindow3D::Text(double x0, double y0, double z0, char *str)
    {
        PlotWindow2D::Text(projx(x0,y0,z0),projy(x0,y0,z0),str);
    }

//
//      PlotWindow3D::Dot
//
//      Print a dot at the physical coordinates with the 
//      specified color (see PlotWindow.cp)
//
//      Created 11/21/95
//
    void PlotWindow3D::Dot(double x0, double y0, double z0, int color)
    {
        PlotWindow2D::Dot(projx(x0,y0,z0), projy(x0,y0,z0), color);
    }

//
//      PlotWindow3D::Triangle
//
//      Plot a triangle with vertices at the given physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow3D::Triangle(double x0, double y0, double z0,
                                double x1, double y1, double z1,
                                double x2, double y2, double z2)
    {
        double vec[3];
        double norm;
        unsigned short shade;

        vec[0] = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
        vec[1] = (x2-x0)*(z1-z0)-(x1-x0)*(z2-z0);
        vec[2] = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
        norm = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
        shade = (unsigned short)(USHRT_MAX*fabs((vec[0]*view[0]+vec[1]*view[1]
                                                 +vec[2]*view[2])/norm/sqrt(vdv)));
        PlotWindow2D::Triangle(projx(x0,y0,z0), projy(x0,y0,z0),
                               projx(x1,y1,z1), projy(x1,y1,z1),
                               projx(x2,y2,z2), projy(x2,y2,z2), shade);
    }

//
//      PlotWindow3D::Polygon
//
//      Plot a polygon with vertices at the given physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow3D::Polygon(int num, double* x, double* y, double* z,
                               unsigned short color)
    {
        double *tx = new double[num];
        double *ty = new double[num];
        for (int i=0; i<num; ++i) {
            tx[i] = projx(x[i],y[i],z[i]);
            ty[i] = projy(x[i],y[i],z[i]);
        }
        PlotWindow2D::Polygon(num, tx, ty, color);
        delete[] tx;
        delete[] ty;
    }

//
//      PlotWindow3D::Circle
//
//      Plot a circle centered at (x0,y0) with radius r in the physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow3D::Circle(double x0, double y0, double z0, double radius)
    {
        PlotWindow2D::Circle(projx(x0,y0,z0), projy(x0,y0,z0), radius);
    }

//
//      PlotWindow3D::FilledCircle
//
//      Plot a circle centered at (x0,y0) with radius r in the physical
//      coordinates
//
//      Created 11/21/95
//
    void PlotWindow3D::FilledCircle(double x0, double y0, double z0, double radius)
    {
        PlotWindow2D::FilledCircle(projx(x0,y0,z0), projy(x0,y0,z0), radius);
    }

}
