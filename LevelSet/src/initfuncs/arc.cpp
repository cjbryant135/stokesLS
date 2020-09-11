/*************************************************
 arc.cpp
 
 $Header: arc.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $
 
 $Log:  arc.cpp,v $
 Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
 *** none ***
 
 Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
 Initial revision
 
 Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
 Initial revision
 
 
*************************************************/

#include "arc.h"
#include "utility.h"
#include <limits.h>
#define A parameter[Height]
#define B parameter[Width]

namespace levelset {
        
    void Arc::SetParams(InputParams *params)
    {
        if (parameter) delete[] parameter;
        parameter = new double[2];
        A = params->GetDoubleParam("Height");
        B = params->GetDoubleParam("xmax");
    }
        
    double Arc::XY(const double s, const double t) const
    {
        double answer;
        if (A == 0.) {
            answer = t;
        } else {
            double c = (B*B-A*A)/2./A;
            answer = sqrt(B*B+c*c)-sqrt(s*s+(t-c)*(t-c));
        }
        return answer;
    }
        
}






