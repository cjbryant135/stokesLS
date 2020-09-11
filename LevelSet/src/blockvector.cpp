#include "blockvector.h"
#include <math.h>
#include <iostream>
#include <iomanip>

namespace blockmatrix {

    BlockVector::BlockVector(const size_t dim) : Dim(dim), vec(NULL),
                                                 Alldim(0)
    {
        vec = new spVector*[dim];
        for (int r=0; r<dim; ++r) vec[r] = NULL;
    }

    void BlockVector::SubVector(const size_t r, spVector& v)
    {
        if (vec[r] != NULL) Alldim -= v.length();
        vec[r] = &v;
        Alldim += v.length();
    }

    BlockVector::BlockVector(const BlockVector& W) 
        : Dim(W.Dim), vec(NULL), Alldim(W.Alldim)
    {
        vec = new spVector*[Dim];
        for (int r=0; r<Dim; ++r) vec[r] = new spVector(*W.vec[r]);
    }

    BlockVector::~BlockVector(void)
    {
        if (vec != NULL)
            delete[] vec;
    }

    double BlockVector::operator()(const size_t i) const 
    {
        int r, j;
        for (r=0, j=i; r<Dim && j>=vec[r]->length(); j-=vec[r++]->length()) ;
        return (*vec[r])(j);
    }

    double& BlockVector::operator()(const size_t i)
    {
        int r, j;
        for (r=0, j=i; r<Dim && j>=vec[r]->length(); j-=vec[r++]->length()) ;
        return (*vec[r])(j);
    }

    void BlockVector::operator=(const BlockVector& W)
    {
        for (int r=0; r<Dim; ++r)
            *vec[r] = *W.vec[r];
    }

    void BlockVector::operator=(const double a)
    {
        for (int r=0; r<Dim; ++r)
            *vec[r] = a;
    }

    void BlockVector::operator+=(const BlockVector& W)
    {
        for (int r=0; r<Dim; ++r)
            *vec[r] += *W.vec[r];
    }
   
    void BlockVector::operator*=(const double a)
    {
        for (int r=0; r<Dim; ++r)
            *vec[r] *= a;
    }

    void BlockVector::operator*=(const double* W)
    {
        for (int r=0, i=0; r<Dim; ++r) {
            *vec[r] *= &(W[i-1]);
            i += vec[r]->length();
        }
    }
   
    void BlockVector::operator-=(const BlockVector& W)
    {
        for (int r=0; r<Dim; ++r)
            *vec[r] -= *W.vec[r];
    }

    void BlockVector::operator-=(const double* v)
    {
        for (int r=0, i=0; r<Dim; ++r) {
            *vec[r] -= &(v[i]);
            i += vec[r]->length();
        }
    }

    void BlockVector::operator/=(const double a)
    {
        for (int r=0; r<Dim; ++r)
            *vec[r] /= a;
    }

    double BlockVector::l1Norm(void) const
    {
        double Sum;
   
        Sum = 0.0;
        for (int r=0; r<Dim; ++r)
            Sum += vec[r]->l1Norm();
   
        return(Sum);
    }

    double BlockVector::l2Norm(void) const
    {
        double Sum;
        double temp;
   
        Sum = 0.0;
        for (int r=0; r<Dim; ++r) {
            temp = vec[r]->l2Norm();
            Sum += temp*temp;
        }
   
        return(sqrt(Sum));
    }

    double BlockVector::linfNorm(void) const
    {
        double mx;
        double tmp;
   
        mx = 0.0;
        for (int r=0; r<Dim; ++r) {
            tmp = vec[r]->linfNorm();
            mx = mx > tmp ? mx : tmp;
        }
   
        return(mx);
    }

    void BlockVector::print(std::ostream& s) const
    {
        for (int r=0; r<Dim; ++r) {
            s << "Block " << r << '\n';
            vec[r]->print(s);
        }
    }

    void BlockVector::matlab_print(const char* n, std::ostream& s) const
    {
        s << n << " = [\n";
        for (int r=0; r<Dim; ++r) {
            for (int i=1; i<=vec[r]->Dim; ++i) {
                s << std::setiosflags(std::ios::scientific) << vec[r]->Cmp[i] << ";\n";
            }
        }
        s << "];\n";
        s.flush();
    }

}
