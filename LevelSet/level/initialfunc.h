#ifndef __INITIALFUNC_H__
#define __INITIALFUNC_H__

#include "defs.h"
#include "inputparams.h"

namespace levelset {

    typedef double (*DDfunc)(double,double);
    typedef float (*FFfunc)(float,float);
    typedef int (*IIfunc)(int,int);
    typedef double (*Dfunc)(double);
    typedef float (*Ffunc)(float);
    typedef int (*Ifunc)(int);

    class InitialFunc {
    public:

    InitialFunc(void) : start(0), stop(1), xy_numsteps(1000), parameter(NULL) {}

    InitialFunc(const int nparams, const double* params) 
        : start(0), stop(1), xy_numsteps(1000), parameter(NULL) {SetParams(nparams, params);} 

    InitialFunc(const double p1) : start(0), stop(1), xy_numsteps(1000), parameter(NULL)
        {SetParams(p1);}

    InitialFunc(const double p1, const double p2) : start(0), stop(1),
            xy_numsteps(1000), parameter(NULL)
        {SetParams(p1, p2);} 

    InitialFunc(const double p1, const double p2, const double p3)
        : start(0), stop(1), xy_numsteps(1000), parameter(NULL) 
        {SetParams(p1, p2, p3);}
   
    InitialFunc(const double p1, const double p2, const double p3,
                const double p4) : start(0), stop(1), xy_numsteps(1000), parameter(NULL)
        {SetParams(p1, p2, p3, p4);}

    InitialFunc(InputParams *params) : start(0), stop(1), xy_numsteps(1000), parameter(NULL)
        {SetParams(params);} 
   
        virtual ~InitialFunc(void) {if (parameter) delete[] parameter;}

        virtual void SetParams(InputParams *params) {} 
   
        virtual double X(const double t) const {return 0;}
        virtual double Y(const double t) const {return 0;}
        virtual double Z(const double t) const {return 0;}
        virtual double XY(const double s, const double t) const;
        virtual double XYZ(const double s, const double t, const double u) const
        {return XY(s,t);}

        double Start(void) const {return start;}
        double Stop(void) const {return stop;}

    protected:
        int     param_num;
        double       start, stop;
        const int xy_numsteps;
        double*      parameter;

        void SetParams(const int nparams, const double* params);
        void SetParams(const double p1);
        void SetParams(const double p1, const double p2);
        void SetParams(const double p1, const double p2, const double p3);
        void SetParams(const double p1, const double p2, const double p3,
                       const double p4);
   
        double Brent(const double bx, const double ds, const double s, const double t) const;
    };

}

#endif //__INITIALFUNC_H__
