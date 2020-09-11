/*************************************************
    mcubes.cpp

    $Header: mcubes.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $

    $Log:       mcubes.cpp,v $
Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
Initial revision


*************************************************/

#include "mcubes.h"
#include "utility.h"
#include <float.h>
#include <sstream>

namespace levelset {

    void Mcubes::SetParams(InputParams *params)
    {
        //char s[255];
        num = params->GetIntParam("Number of cubes");
        if (parameter) delete[] parameter;
        parameter = new double[4*num];
        for (int i=0; i<num; ++i) {
            {
                std::ostringstream str;
                str << "Radius of Cube " << i+1 << std::ends;
                parameter[4*i] = params->GetDoubleParam(str.str().c_str());
            }
            {
                std::ostringstream str;
                str << "X Center of Cube " << i+1 << std::ends;
                parameter[4*i+1] = params->GetDoubleParam(str.str().c_str());
            }
            {
                std::ostringstream str;
                str << "Y Center of Cube " << i+1 << std::ends;
                parameter[4*i+2] = params->GetDoubleParam(str.str().c_str());
            }
            {
                std::ostringstream str;
                str << "Z Center of Cube " << i+1 << std::ends;
                parameter[4*i+3] = params->GetDoubleParam(str.str().c_str());
            }
        }
                
    }

    double Mcubes::XYZ(const double s, const double t, const double u) const
    {
        double answer = -DBL_MAX;
        for (int i=0; i<num; ++i)
            answer = max(answer,parameter[4*i]-max(fabs(s-parameter[4*i+1]),
                                                   fabs(t-parameter[4*i+2]),fabs(u-parameter[4*i+3])));
   
        return answer;
    }

}







