#include <iostream>
#include "blockqmatrix.h"

namespace blockmatrix {
        
    BlockQMatrix::BlockQMatrix(const size_t rc) : Dim(rc), Alldim(0),
                                                  mat(NULL), InvDiagEl(NULL)
    {
        mat = new spMat*[Dim];
        for (int r=0; r<Dim; ++r) {
            mat[r] = new spMat[Dim]; 
            for (int c=0; c<Dim; ++c)
                mat[r][c].spq = NULL;
        }
    }

    BlockQMatrix::~BlockQMatrix(void)
    {
        Clear();
    }

    void BlockQMatrix::Clear(const size_t r, const size_t c)
    {
        if (r == UINT_MAX || c == UINT_MAX) {
            for (int r=0; r<Dim; ++r)
                for (int c=0; c<Dim; ++c)
                    if (mat[r][c].spq != NULL) {
                        if (r == c)
                            mat[r][c].spq->Clear();
                        else
                            mat[r][c].sp->Clear();
                    }
        }
        else {
            if (r == c)
                mat[r][c].spq->Clear();
            else
                mat[r][c].sp->Clear();
        }
        if (InvDiagEl != NULL) {
            delete[] InvDiagEl;
            InvDiagEl = NULL;
        }
    }

    void BlockQMatrix::SubMatrix(const size_t rc, spQMatrix& m)
    {
        if (mat[rc][rc].spq != NULL) Alldim -= mat[rc][rc].spq->size();
        mat[rc][rc].spq = &m;
        Alldim += mat[rc][rc].spq->size();
    }

    void BlockQMatrix::SubMatrix(const size_t r, const size_t c, spMatrix& m) 
    {
        mat[r][c].sp = &m;
    }

    void BlockQMatrix::FindSubMatrix(const size_t r, const size_t c, 
                                     size_t& i, size_t& j,
                                     size_t& rr, size_t& cc) const
    {
        i = 0;
        if (r >= mat[0][0].spq->size()) 
            for (rr = r-mat[0][0].spq->size(), i=1; i < Dim-1 && rr >= mat[0][i].sp->cols();
                 rr -= mat[0][i++].sp->cols()) ;
        else 
            rr = r;
        j = 0;
        if (c >= mat[0][0].spq->size())
            for (cc = c-mat[0][0].spq->size(), j=1; j < Dim-1 && cc >= mat[0][j].sp->cols();
                 cc -= mat[0][j++].sp->cols()) ;
        else
            cc = c;
    }

    void BlockQMatrix::SortEl(void)
    {
        for (int i=0; i<Dim; ++i)
            mat[i][i].spq->SortEl();
    }

    void BlockQMatrix::AllocInvDiagEl(void)
    {
        if (InvDiagEl == NULL) 
            InvDiagEl = new double[Alldim];
        for (int i=0, j=0; i<Dim; ++i) {
            mat[i][i].spq->AllocInvDiagEl();
            for (int k=0; k<mat[i][i].spq->Dim; ++k, ++j)
                InvDiagEl[j] = mat[i][i].spq->InvDiagEl[k+1];
        }
    }

    double BlockQMatrix::operator()(const size_t r, const size_t c) const
    {
        return Value(r,c);
    }

    double BlockQMatrix::Value(const size_t r, const size_t c) const
    {
        size_t i, j, rr, cc;
        
        FindSubMatrix(r, c, i, j, rr, cc);
        return(i == j ? mat[i][i].spq->Value(rr,cc) 
               : mat[i][j].sp->Value(rr,cc));
    }

    double& BlockQMatrix::operator()(const size_t r, const size_t c)
    {
        size_t i, j, rr, cc;
        FindSubMatrix(r, c, i, j, rr, cc);
        if (i == j)
            return (*mat[i][j].spq)(rr,cc);
        else
            return (*mat[i][j].sp)(rr,cc);
    }

    void BlockQMatrix::Zero(const size_t r, const size_t c)
    {
        size_t i, j, rr, cc;
        FindSubMatrix(r, c, i, j, rr, cc);
        if (i == j)
            return mat[i][j].spq->Zero(rr,cc);
        else
            return mat[i][j].sp->Zero(rr,cc);
    }

    double* BlockQMatrix::ColSums(double* s)
    {
        double *sum, *temp;

        if (s == NULL) {
            sum = new double[Alldim];
        } else
            sum = s;
        for (int r=0; r<Alldim; ++r) sum[r] = 0.;

        temp = new double[Alldim];
        for (int r=0; r<Dim; ++r) 
            for (int c=0, i=0; c<Dim; ++c) 
                if (r == c) {
                    mat[r][c].spq->ColSums(temp);
                    for (int k=0; k<mat[r][c].spq->cols(); ++k)
                        sum[i+k] += temp[k];
                    i += mat[r][c].spq->cols();
                } else {
                    mat[r][c].sp->ColSums(temp);
                    for (int k=0; k<mat[r][c].sp->cols(); ++k)
                        sum[i+k] += temp[k];
                    i += mat[r][c].sp->cols();
                }
   
        delete[] temp;
        return sum;
    }     

    BlockQMatrix& BlockQMatrix::operator=(const BlockQMatrix& A)
    {
        for (int r=0; r<Dim; ++r)
            for (int c=0; c<Dim; ++c)
                if (r==c)
                    *mat[r][c].spq = *A.mat[r][c].spq;
                else
                    *mat[r][c].sp = *A.mat[r][c].sp;
        return *this;
    }

    void BlockQMatrix::print(std::ostream& s) const
    {
        for (int r=0; r<Dim; ++r)
            for (int c=0; c<Dim; ++c) {
                s << "Matrix (" << r << ',' << c << ")\n";
                if (r==c) 
                    mat[r][c].spq->print(s);
                else
                    mat[r][c].sp->print(s);
            }
    }

    void BlockQMatrix::matlab_print(const char* n, std::ostream& s) const
    {
        s << n << " = zeros(" << Alldim << ',' << Alldim << ");\n";
        for (int r=0, k=0; r<Dim; ++r) {
            for (int c=0, l=0; c<Dim; ++c) {
                if (r==c) {
                    for (int i=1; i<=mat[r][c].spq->Dim; ++i) {
                        for (int j=0; j<mat[r][c].spq->Len[i]; ++j)
                            if (mat[r][c].spq->El[i][j].Pos)
                                s << n << "(" << k+i << ',' 
                                  << l+mat[r][c].spq->El[i][j].Pos << ") = " 
                                  << mat[r][c].spq->El[i][j].Val << ";\n";
                    }
                    l += mat[r][c].spq->cols();
                }
                else {
                    for (int i=1; i<=mat[r][c].sp->RowDim; ++i) {
                        for (int j=0; j<mat[r][c].sp->Len[i]; ++j)
                            if (mat[r][c].sp->El[i][j].Pos)
                                s << n << "(" << k+i << ',' 
                                  << l+mat[r][c].sp->El[i][j].Pos << ") = " 
                                  << mat[r][c].sp->El[i][j].Val << ";\n";
                    }
                    l += mat[r][c].sp->cols();
                }
            }
            k += mat[r][r].spq->rows();
        }
        s.flush();
    }

}
        

