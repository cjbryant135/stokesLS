/*************************************************
    spvector.h

    $Header: spvector.h,v 1.1 99/09/20 11:36:12 chopp Exp $

    $Log:       spvector.h,v $
    * Revision 1.1  99/09/20  11:36:12  11:36:12  chopp (David Chopp)
    * Initial revision
    * 
    * Revision 1.3  99/04/05  13:32:35  13:32:35  chopp (David Chopp)
    * *** none ***
    * 
    * Revision 1.2  99/03/02  11:53:26  11:53:26  chopp (David Chopp)
    * *** none ***
    * 
    * Revision 1.1  99/02/26  14:26:40  14:26:40  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __SPVECTOR_H__
#define __SPVECTOR_H__

#include <stdlib.h>
#include <iostream>

namespace blockmatrix {
        
    class spVector 
    {
    public:

        size_t Dim;
        double *Cmp;

        spVector(const size_t dim = 1);
        spVector(const spVector& W);
        ~spVector(void);

        void reshape(const size_t i);
        void resize(const size_t i);
        size_t length(void) const {return Dim;}

        double l1Norm(void) const;
        double l2Norm(void) const;
        double linfNorm(void) const;

        inline double operator()(const size_t i) const {return Cmp[i+1];}
        inline double& operator()(const size_t i) {return Cmp[i+1];}

        void operator=(const spVector& W);
        void operator=(const double a);
        void operator+=(const spVector& W);
        void operator*=(const double a);
        void operator*=(const double* d);
        void operator-=(const spVector& W);
        void operator-=(const double* d);
        void operator/=(const double a);

        void print(std::ostream& s = std::cout) const;
        void matlab_print(const char* n, std::ostream& s = std::cout) const;
   
    };

    double operator*(const spVector& V, const spVector& W);

}
#endif /* __SPVECTOR_H__ */

