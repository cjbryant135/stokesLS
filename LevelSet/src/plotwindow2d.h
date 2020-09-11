
#ifdef __SC__
#include <Quickdraw.h>
#endif
#include "plotwindow.h"

#ifndef H_PLOTWINDOW2D_CLASS
#define H_PLOTWINDOW2D_CLASS

namespace levelset {
        
    class PlotWindow2D : public PlotWindow {
    protected:
   
        double  pxmin, pxmax, pymin, pymax; // Physical dimensions of the plot region
        
        inline int ispot(double x) {return int((x-pxmin)/(pxmax-pxmin)*10000.00);};
        inline int jspot(double y) {return int((y-pymin)/(pymax-pymin)*10000.00);};

    public:

        // PlotWindow2D(const int w = 800, const int h = 800) : PlotWindow(w,h) {;}
    PlotWindow2D(bool yes_window = true, bool pae = true, const int w = 800, const int h = 800)
        : PlotWindow(yes_window, pae, w, h) {;}
        virtual ~PlotWindow2D(void);
        
        virtual void SetRange(double xmin, double xmax, double ymin, double ymax);
        virtual void Line(double x0, double y0, double x1, double y1);
        virtual void Arrow(double x0, double y0, double x1, double y1);
        virtual void Text(double x0, double y0, const char *str);
//   void Text(WinPos wp, char* str);
        virtual void Dot(double x0, double y0, int color = 0);
        virtual void Triangle(double x0, double y0, double x1, double y1, double x2,
                              double y2, unsigned short color);
        virtual void Polygon(int num, double* x, double* y, unsigned short color);
        virtual void Circle(double x0, double y0, double radius);
        virtual void FilledCircle(double x0, double y0, double radius);
   
        virtual void Clear(double xmin, double xmax, double ymin, double ymax);
        virtual void Clear(const bool pause = false) {PlotWindow::Clear(pause);}
    };

}
#endif
