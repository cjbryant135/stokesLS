/*************************************************
    mballs.cpp

    $Header: mballs.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $

    $Log:       mballs.cpp,v $
Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
Initial revision


*************************************************/

#include "mballs.h"
#include "utility.h"
#include <float.h>
#include <sstream>

namespace levelset {

    void Mballs::SetParams(InputParams *params)
    {
        num = params->GetIntParam("Number of Balls",1,"");
        if (parameter) delete[] parameter;
        parameter = new double[4*num];
        for (int i=0; i<num; ++i) {
            {
                std::ostringstream str;
                str << "Radius of Ball " << i+1 << std::ends;
                parameter[4*i] = params->GetDoubleParam(str.str().c_str(),1.,"");
            }
            {
                std::ostringstream str;
                str << "X Center of Ball " << i+1 << std::ends;
                parameter[4*i+1] = params->GetDoubleParam(str.str().c_str(),0.,"");
            }
            {
                std::ostringstream str;
                str << "Y Center of Ball " << i+1 << std::ends;
                parameter[4*i+2] = params->GetDoubleParam(str.str().c_str(),0.,"");
            }
            {
                std::ostringstream str;
                str << "Z Center of Ball " << i+1 << std::ends;
                parameter[4*i+3] = params->GetDoubleParam(str.str().c_str(),0.,"");
            }
        }
                
    }

    double Mballs::XYZ(const double s, const double t, const double u) const
    {
        double answer = -DBL_MAX;
        for (int i=0; i<num; ++i)
            answer = max(answer,parameter[4*i]-sqrt(sqr(s-parameter[4*i+1])
                                                    +sqr(t-parameter[4*i+2])+sqr(u-parameter[4*i+3])));
   
        return answer;
    }


}






