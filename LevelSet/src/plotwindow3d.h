/*************************************************
    plotwindow3d.h

    $Header: plotwindow3d.h,v 1.1 99/02/04 14:38:00 chopp Exp $

    $Log:       plotwindow3d.h,v $
    * Revision 1.1  99/02/04  14:38:00  14:38:00  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/


#ifdef __SC__
#include <Quickdraw.h>
#endif
#include "plotwindow2d.h"

#ifndef H_PLOTWINDOW3D_CLASS
#define H_PLOTWINDOW3D_CLASS

namespace levelset {
        
    class PlotWindow3D : public PlotWindow2D {

// Physical dimensions of the plot region
        double limits[3][2];
        double view[3];
        double a_[2], b_[3];
        double vdv;

        inline double projx(double x, double y, double z)
        {
            return((a_[0]*x+a_[1]*y)/(vdv-x*view[0]-y*view[1]-z*view[2]));
        }
   
        inline double projy(double x, double y, double z)
        {
            return((b_[0]*x+b_[1]*y+b_[2]*z)/(vdv-x*view[0]-y*view[1]-z*view[2]));
        }

    public:

        //PlotWindow3D(const int w = 800, const int h = 800) : PlotWindow2D(w,h) {;}
    PlotWindow3D(bool yes_window = true, bool pae = true, const int w = 800, const int h = 800)
        : PlotWindow2D(yes_window, pae, w, h) 
        {view[0] = -1000.; view[1] = -2000; view[2] = 1500;}
        virtual ~PlotWindow3D(void);
        
        virtual void SetRange(double xmin, double xmax, double ymin, double ymax,
                              double zmin, double zmax);
        virtual void GetRange(double& xmin, double& xmax, double& ymin, double& ymax,
                              double& zmin, double& zmax) const;
        virtual void GetActualRange(double& xmin, double& xmax, double& ymin, double& ymax,
                                    double& zmin, double& zmax) const;
        virtual void SetView(double x, double y, double z);
        virtual void GetView(double& x, double& y, double& z) const;
        virtual void Axes(void);
        virtual void Line(double x0, double y0, double z0, double x1, double y1, double z1);
        virtual void Arrow(double x0, double y0, double z0,
                           double x1, double y1, double z1);
        virtual void Text(double x0, double y0, double z0, char *str);
        virtual void Dot(double x0, double y0, double z0, int color = 0);
        virtual void Triangle(double x0, double y0, double z0,
                              double x1, double y1, double z1,
                              double x2, double y2, double z2);
        virtual void Polygon(int num, double* x, double* y, double* z,
                             unsigned short color);
        virtual void Circle(double x0, double y0, double z0, double radius);
        virtual void FilledCircle(double x0, double y0, double z0, double radius);

    private:

        void Recalc(void);
    };

}
#endif
