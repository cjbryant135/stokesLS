/*************************************************
    mcylinders.cpp

    $Header: mcylinders.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $

    $Log:       mcylinders.cpp,v $
Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
Initial revision


*************************************************/

#include "mcylinders.h"
#include "utility.h"
#include <float.h>
#include <sstream>

namespace levelset {
        
    void Mcylinders::SetParams(InputParams *params)
    {
        num = params->GetIntParam("Number of Cylinders",1," ");
        if (parameter) delete[] parameter;
        parameter = new double[dpc*num];
        for (int i=0; i<num; ++i) {
            {
                std::ostringstream str;
                str << "Radius of Cylinder " << i+1 << std::ends;
                parameter[dpc*i] = params->GetDoubleParam(str.str().c_str(),1.,"");
            }
            {
                std::ostringstream str;
                str << "X1 End of Cylinder " << i+1 << std::ends;
                parameter[dpc*i+1] = params->GetDoubleParam(str.str().c_str(),0.,"");
            }
            {
                std::ostringstream str;
                str << "Y1 End of Cylinder " << i+1 << std::ends;
                parameter[dpc*i+2] = params->GetDoubleParam(str.str().c_str(),0.,"");
            }
            {
                std::ostringstream str;
                str << "Z1 End of Cylinder " << i+1 << std::ends;
                parameter[dpc*i+3] = params->GetDoubleParam(str.str().c_str(),0.,"");
            }
            {
                std::ostringstream str;
                str << "X2 End of Cylinder " << i+1 << std::ends;
                parameter[dpc*i+4] = params->GetDoubleParam(str.str().c_str(),
                                                            parameter[dpc*i+1],"");
            }
            {
                std::ostringstream str;
                str << "Y2 End of Cylinder " << i+1 << std::ends;
                parameter[dpc*i+5] = params->GetDoubleParam(str.str().c_str(),
                                                            parameter[dpc*i+2],"");
            }
            {
                std::ostringstream str;
                str << "Z2 End of Cylinder " << i+1 << std::ends;
                parameter[dpc*i+6] = params->GetDoubleParam(str.str().c_str(),
                                                            parameter[dpc*i+3],"");
            }
        }
                
    }

#define R parameter[dpc*i]
#define X1 parameter[dpc*i+1]
#define Y1 parameter[dpc*i+2]
#define Z1 parameter[dpc*i+3]
#define X2 parameter[dpc*i+4]
#define Y2 parameter[dpc*i+5]
#define Z2 parameter[dpc*i+6]

    double area(const double x1, const double y1, const double z1,
                const double x2, const double y2, const double z2)
    {
        double x = y1*z2-y2*z1;
        double y = x2*z1-x1*z2;
        double z = x1*y2-x2*y1;
        return sqrt(x*x+y*y+z*z);
    }

    double Mcylinders::XYZ(const double s, const double t, const double u) const
    {
        double answer = -DBL_MAX;
        for (int i=0; i<num; ++i) {
            double h = area(X2-X1,Y2-Y1,Z2-Z1,s-X1,t-Y1,u-Z1)
                /sqrt(sqr(X2-X1)+sqr(Y2-Y1)+sqr(Z2-Z1)); 
            answer = max(answer,R-h);
        }
        return answer;
    }

}






