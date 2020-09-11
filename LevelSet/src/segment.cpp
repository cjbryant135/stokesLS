#include "segment.h"

namespace levelset {
        
    Segment::Segment(void)
    {
        nabor[0] = -1;
        nabor[1] = -1;
        sides = SNone;
        data = NULL;
    }

    Segment::Segment(const double x1, const double y1, const double x2,
                     const double y2, const int s, const int dl)
    {
        x[0] = x1;
        x[1] = x2;
        y[0] = y1;
        y[1] = y2;
        nabor[0] = -1;
        nabor[1] = -1;
        sides = s;
        data = new double[dl];
    }

    Segment::Segment(const Segment& s)
    {
        x[0] = s.x[0];
        x[1] = s.x[1];
        y[0] = s.y[0];
        y[1] = s.y[1];
        datalen = s.datalen;
        data = new double[datalen];
        for (int i=0; i<datalen; ++i) data[i] = s.data[i];
        nabor[0] = s.nabor[0];
        nabor[1] = s.nabor[1];
        sides = s.sides;
    }

    Segment::~Segment(void)
    {
        if (data != NULL) delete[] data;
    }

    void Segment::SetXY(const double x1, const double y1, const double x2,
                        const double y2, const int s, const int dl)
    {
        x[0] = x1;
        x[1] = x2;
        y[0] = y1;
        y[1] = y2;
        sides = s;
        if (data==NULL) {
            datalen = dl;
            data = new double[datalen];
        }
    }

    Segment& Segment::operator=(const Segment& a)
    {
        x[0] = a.x[0];
        x[1] = a.x[1];
        y[0] = a.y[0];
        y[1] = a.y[1];
        if (a.data != NULL) {
            if (data != NULL) delete[] data;
            datalen = a.datalen;
            data = new double[datalen];
        }
        else {
            if (data != NULL) delete[] data;
            datalen = 0;
            data = NULL;
        }
        for (int ii=0; ii<datalen; ++ii) data[ii] = a.data[ii];
        nabor[0] = a.nabor[0];
        nabor[1] = a.nabor[1];
        sides = a.sides;
        return *this;
    }


#ifndef NO_GRAPHICS
    void Segment::Plot(PlotWindow2D& port, const int vel, const char showvel,
                       const double scale)
    {
        port.Line(x[0], y[0], x[1], y[1]);
        if (showvel) {
            double mx, my;
            MidPoint(mx, my);
            double nx, ny;
            Normal(nx, ny);
            port.Arrow(mx, my, mx+nx*scale*data[vel],
                       my+ny*scale*data[vel]);
        }
    }
#endif
    
}









