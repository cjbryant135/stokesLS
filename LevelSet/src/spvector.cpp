/*************************************************
    spvector.cpp

    $Header: spvector.cpp,v 1.1 99/09/20 11:35:50 chopp Exp $

    $Log:       spvector.cpp,v $
Revision 1.1  99/09/20  11:35:50  11:35:50  chopp (David Chopp)
Initial revision

Revision 1.3  99/04/05  13:32:20  13:32:20  chopp (David Chopp)
*** none ***

Revision 1.2  99/03/02  11:53:23  11:53:23  chopp (David Chopp)
*** none ***

Revision 1.1  99/02/26  14:02:39  14:02:39  chopp (David Chopp)
Initial revision

Revision 1.1  99/02/26  14:01:51  14:01:51  chopp (David Chopp)
Initial revision

*************************************************/

#include "spvector.h"
#include <math.h>
#include <iostream>
#include <iomanip>

namespace blockmatrix {
        
    spVector::spVector(const size_t dim)
    {
        Dim = dim;
        Cmp = new double[Dim+1];
    }

    spVector::spVector(const spVector& W)
    {
        Dim = W.length();
        Cmp = new double[Dim+1];
        for (int d=1; d <= Dim; ++d)
            Cmp[d] = W.Cmp[d];
    }

    spVector::~spVector(void)
    {
        if (Cmp != NULL)
            delete[] Cmp;
    }

    void spVector::reshape(const size_t dim)
    {
        if (Cmp != NULL)
            delete[] Cmp;
        Dim = dim;
        Cmp = new double[Dim+1];
        for (int i=0; i<=Dim; ++i) Cmp[i] = 0.;
    }

    void spVector::resize(const size_t dim)
    {
        double* tmp = new double[dim+1];
        int mindim = dim < Dim ? dim : Dim;
        if (Cmp != NULL) {
            for (int i=1; i<=mindim; ++i)
                tmp[i] = Cmp[i];
            delete[] Cmp;
        }
        Dim = dim;
        Cmp = tmp;
    }

    void spVector::operator=(const spVector& W)
    {
        size_t dim = W.length();
        if (Dim != dim) reshape(dim);
        for (int d=1; d<=Dim; ++d)
            Cmp[d] = W.Cmp[d];
    }

    void spVector::operator=(const double a)
    {
        for (int d=1; d<=Dim; ++d)
            Cmp[d] = a;
    }

    void spVector::operator+=(const spVector& W)
    {
        for (int d=1; d<=Dim; ++d)
            Cmp[d] += W.Cmp[d];
    }
   
    void spVector::operator*=(const double a)
    {
        for (int d=1; d<=Dim; ++d)
            Cmp[d] *= a;
    }

    void spVector::operator*=(const double* W)
    {
        for (int d=1; d<=Dim; ++d)
            Cmp[d] *= W[d];
    }
   
    void spVector::operator-=(const spVector& W)
    {
        for (int d=1; d<=Dim; ++d)
            Cmp[d] -= W.Cmp[d];
    }

    void spVector::operator-=(const double* v)
    {
        for (int d=1; d<=Dim; ++d) 
            Cmp[d] -= v[d-1];
    }

    void spVector::operator/=(const double a)
    {
        for (int d=1; d<=Dim; ++d)
            Cmp[d] /= a;
    }

    double spVector::l1Norm(void) const
    {
        double Sum;
   
        Sum = 0.0;
        for (int i=1; i<=Dim; ++i)
            Sum += fabs(Cmp[i]);
   
        return(Sum);
    }

    double spVector::l2Norm(void) const
    {
        double Sum;
   
        Sum = 0.0;
        for (int i=1; i<=Dim; ++i)
            Sum += Cmp[i]*Cmp[i];
   
        return(sqrt(Sum));
    }

    double spVector::linfNorm(void) const
    {
        double mx;
   
        mx = 0.0;
        for (int i=1; i<=Dim; ++i)
            mx = mx > fabs(Cmp[i]) ? mx : fabs(Cmp[i]);
   
        return(mx);
    }

    void spVector::print(std::ostream& s) const
    {
        int lcount = 5;
        for (int i=1; i<=Dim; ++i) {
            s << std::setiosflags(std::ios::scientific) << Cmp[i] << '\t';
            if (i%lcount == 0)
                s << '\n';
        }
    }

    void spVector::matlab_print(const char* n, std::ostream& s) const
    {
        s << n << " = [\n";
        for (int i=1; i<=Dim; ++i) {
            s << std::setiosflags(std::ios::scientific) << Cmp[i] << ";\n";
        }
        s << "];\n";
        s.flush();
    }

}
