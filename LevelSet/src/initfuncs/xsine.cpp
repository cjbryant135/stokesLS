/*************************************************
    xsine.cpp

    $Header: xsine.cpp,v 1.1 2000/05/31 11:07:56 chopp Exp $

    $Log:       xsine.cpp,v $
Revision 1.1  2000/05/31  11:07:56  11:07:56  chopp (David Chopp)
Initial revision

*************************************************/

#include "xsine.h"
#include "utility.h"
#define A parameter[Amplitude]
#define B parameter[Height]
#define W parameter[Frequency]
#define D parameter[Offset]

namespace levelset {

    XSinusoidal::XSinusoidal(const  double amp, const double ht, const double freq, const double off)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[4];
        A = amp;
        B = ht;
        W = freq;
        D = off;
    }

    void XSinusoidal::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[4];
        A = param->GetDoubleParam("Amplitude",0.,"Amplitude of wave");
        B = param->GetDoubleParam("Height",0.,"Mean Y position of wave");
        W = param->GetDoubleParam("Frequency",0.,"Frequency of wave");
        D = param->GetDoubleParam("Offset",0.,"Offset of wave");
    }

    double XSinusoidal::X(const double t) const
    {
        return t;
    }

    double XSinusoidal::Y(const double t) const
    {
        return A*sin(W*t-D)+B;
    }

    double XSinusoidal::XY(const double s, const double t) const
    {
#if 0
        std::cout << "A = " << A << '\n';
        std::cout << "W = " << W << '\n';
        std::cout << "D = " << D << '\n';
        std::cout << "B = " << B << '\n';
#endif
        int n = (int)((W*s-D)/M_PI+0.5);
        double xm;
        double xl = ((n-0.5)*M_PI+D)/W;
        double xr = ((n+0.5)*M_PI+D)/W;
        double Fl = 2*(xl-s)-2*(t-A*sin(W*xl-D)-B)*A*cos(W*xl-D)*W;
        double Fr = 2*(xr-s)-2*(t-A*sin(W*xr-D)-B)*A*cos(W*xr-D)*W;
        while (xr-xl > 1.0e-10) {
            xm = (xl+xr)/2.;
            double Fm = 2*(xm-s)-2*(t-A*sin(W*xm-D)-B)*A*cos(W*xm-D)*W;
            if (Fl*Fm > 0) {
                xl = xm;
                Fl = Fm;
            } else {
                xr = xm;
                Fr = Fm;
            }
        }
#if 0
        std::cout << "(x0,y0) = (" << s << ", " << t << ")\n";
        std::cout << "(x1,y1) = (" << xm << ", " << A*sin(W*xm-D)+B << ")\n";
        std::cout << "A*sin(W*s-D)+B = " << A*sin(W*s-D)+B << "\n";
#endif
        return sqrt((s-xm)*(s-xm)+(t-A*sin(W*xm-D)-B)*(t-A*sin(W*xm-D)-B))*(t > A*sin(W*s-D)+B ? 1 : -1);
    }

}
