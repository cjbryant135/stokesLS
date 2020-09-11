/*************************************************
    mdisks.cpp

    $Header: mdisks.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $

    $Log:       mdisks.cpp,v $
Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
Initial revision


*************************************************/

#include "mdisks.h"
#include "utility.h"
#include <float.h>
#include <sstream>

namespace levelset {
        
    void Mdisks::SetParams(InputParams *params)
    {
        //char s[255];
        num = params->GetIntParam("Number of Circles");
        if (parameter) delete[] parameter;
        parameter = new double[3*num];
        for (int i=0; i<num; ++i) {
            {
                std::ostringstream str;
                str << "Radius of Circle " << i+1 << std::ends;
                parameter[3*i] = params->GetDoubleParam(str.str().c_str());
            }
            {
                std::ostringstream str;
                str << "X Center of Circle " << i+1 << std::ends;
                parameter[3*i+1] = params->GetDoubleParam(str.str().c_str());
            }
            {
                std::ostringstream str;
                str << "Y Center of Circle " << i+1 << std::ends;
                parameter[3*i+2] = params->GetDoubleParam(str.str().c_str());
            }
        }
                
    }

    void Mdisks::SetParams(const int n, const double *xc, const double *yc, 
                           const double *r)
    {
        num = n;
        if (parameter) delete[] parameter;
        parameter = new double[3*num];
        for (int i=0; i<num; ++i) {
            parameter[3*i] = r[i];
            parameter[3*i+1] = xc[i];
            parameter[3*i+2] = yc[i];
        }
    }



    double Mdisks::XY(const double s, const double t) const
    {
        double answer = -DBL_MAX;
        for (int i=0; i<num; ++i)
            answer = max(answer,parameter[3*i]-sqrt(sqr(s-parameter[3*i+1])
                                                    +sqr(t-parameter[3*i+2])));
   
        return answer;
    }


}






