/*************************************************
    sine.cpp

    $Header: sine.cpp,v 1.1 2000/05/31 11:07:56 chopp Exp $

    $Log:       sine.cpp,v $
Revision 1.1  2000/05/31  11:07:56  11:07:56  chopp (David Chopp)
Initial revision

*************************************************/

#include "sine.h"
#include "utility.h"
#define NX parameter[XNormal]
#define NY parameter[YNormal]
#define B parameter[Hyperplane]
#define A parameter[Amp]
#define F parameter[Freq]
#define O parameter[Offset]

namespace levelset {
        
    Sine::Sine(const  double xn, const double yn, const double b, 
               const double a, const double f, const double off)
        : InitialFunc()
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[6];
        NX = xn;
        NY = yn;
        B = b;
        A = a;
        F = f;
        O = off;
    }

    void Sine::SetParams(InputParams *param)
    {
        start = 0;
        stop = 2.*M_PI;
        if (parameter) delete[] parameter;
        parameter = new double[5];
        NX = param->GetDoubleParam("X-Normal",0.,"X-component of line normal");
        NY = param->GetDoubleParam("Y-Normal",1.,"Y-component of line normal");
        B = param->GetDoubleParam("Hyperplane",0., 
                                  "Line formed by (X-normal,Y-normal).(x,y)=Hyperplane");
        A = param->GetDoubleParam("Perturbation Amplitude",0.1,
                                  "Amplitude of Sine curve perturbation");
        F = param->GetDoubleParam("Perturbation Frequency",1,
                                  "Frequency of Sine curve perturbation");
        O = param->GetDoubleParam("Perturbation Offset",0.,
                                  "Offset of Sine curve perturbation");
    }

    double Sine::X(const double t) const
    {
        return t;
    }

    double Sine::Y(const double t) const
    {
        return 0.;
    }

    double Sine::XY(const double s, const double t) const
    {
        double d = sqrt(NX*NX+NY*NY);
        double ss = (NY*s-NX*t)/d;
        double tt = (NX*s+NY*t)/d-B;
        return tt-A*sin(F*ss-O);
    }

}
