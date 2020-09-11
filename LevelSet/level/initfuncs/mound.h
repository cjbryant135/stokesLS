#ifndef __MOUND_H__
#define __MOUND_H__

#include <math.h>
#include "initialfunc3d.h"

namespace levelset {
        
    class Mound : public InitialFunc3D
    {

        double **x, **y;
        int numc, *num;
        int imax, jmax;
        double xmax, ymax;
        
        void CleanUp(void);
   
    public:

    Mound(InputParams *params) : InitialFunc3D(), x(NULL), y(NULL), numc(0), num(NULL)
        {SetParams(params);} 
        ~Mound(void) {CleanUp();}

        virtual double XYZ(const double s, const double t, const double u) const;

    protected:
        virtual void SetParams(InputParams *params);
   
        void ReadData(const char* fname);
    };

}
#endif
