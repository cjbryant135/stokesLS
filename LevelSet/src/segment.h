#ifndef __SEGMENT_H__
#define __SEGMENT_H__

#include <math.h>
#ifndef NO_GRAPHICS
#include "plotwindow2d.h"
#endif

namespace levelset {
        
    enum SegNabor {SLeft = 0, SRight = 1};
    enum SegData {SVelocity = 0, SLayer, SDFlux, SRFlux, STemp, SLast};
    enum SegSide {SNone = 0, SFromLeft = 1, SFromDown = 2, SFromRight = 4,
                  SFromUp = 8, SToLeft = 16, SToDown = 32, SToRight = 64,
                  SToUp = 128};

    class Segment {

        double    x[2], y[2];
        double    *data;
        int       nabor[2];
        int       sides;
        int       datalen;

    public:
   
        Segment(void);
        Segment(const double x1, const double y1, const double x2, const double y2,
                const int s, const int dl);
        Segment(const Segment& s);
        ~Segment(void);

        // Access routines

        void SetXY(const double x1, const double y1, const double x2,
                   const double y2, const int s, const int dl);

        inline double& X1(void) {return x[0];}
        inline double X1(void) const {return x[0];}
        inline double& Y1(void) {return y[0];}
        inline double Y1(void) const {return y[0];}
        inline double& X2(void) {return x[1];}
        inline double X2(void) const {return x[1];}
        inline double& Y2(void) {return y[1];}
        inline double Y2(void) const {return y[1];}

        inline double& Data(const int i) {return data[i];}
        inline double Data(const int i) const {return data[i];}

        inline int& Nabor(const SegNabor s) {return nabor[s];}
        inline int Nabor(const SegNabor s) const {return nabor[s];}

        inline char Left(void) const {return sides & (SFromLeft | SToLeft);} 
        inline char Right(void) const {return sides & (SFromRight | SToRight);} 
        inline char Down(void) const {return sides & (SFromDown | SToDown);} 
        inline char Up(void) const {return sides & (SFromUp | SToUp);} 
        inline char FromLeft(void) const {return sides & SFromLeft;} 
        inline char FromRight(void) const {return sides & SFromRight;} 
        inline char FromDown(void) const {return sides & SFromDown;} 
        inline char FromUp(void) const {return sides & SFromUp;} 
        inline char ToLeft(void) const {return sides & SToLeft;} 
        inline char ToRight(void) const {return sides & SToRight;} 
        inline char ToDown(void) const {return sides & SToDown;} 
        inline char ToUp(void) const {return sides & SToUp;} 

        // Computational routines

        inline void MidPoint(double& mx, double& my) const 
        {mx = (x[0]+x[1])/2.; my = (y[0]+y[1])/2.;}
        inline void Normal(double& nx, double& ny) const
        {nx = (y[1]-y[0])/Length(); ny = (x[0]-x[1])/Length();} 
        inline double Length(void) const 
        {return sqrt((x[0]-x[1])*(x[0]-x[1])+(y[0]-y[1])*(y[0]-y[1]));}

        // Debugging routines

        inline char Verify(const int sindex = SVelocity) const
        {return finite(data[sindex]);} 

        // Overloaded operators

        Segment& operator=(const Segment& a);
        inline double operator[](const int i) const {return data[i];}
        inline double& operator[](const int i) {return data[i];}

        // I/O stuff

#ifndef NO_GRAPHICS
        void Plot(PlotWindow2D& port, const int vel, const char showvel = 0x00, const double scale = 1);
#endif

        // Friend declarations

        friend class UniformMesh2D;
        friend class Cake;
        friend class Interface;
#ifdef USE_AVS
        friend std::ostream& AVS_Out(std::ostream& s, Interface* m, const char* comment);
#endif
    };

}
#endif // __SEGMENT_H__



