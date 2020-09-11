#include "spmatrix.h"
#include <iostream>

#define DELETE(A) {delete A;A=NULL;}
#define DELARR(A) {delete[] A; A=NULL;}

namespace blockmatrix {
        
    spMatrix::spMatrix(const size_t r, const size_t c) : Len(NULL), El(NULL)
    {
        RowDim = r;
        ClmDim = c;
        Len = new size_t[RowDim+1];
        El = new ElType*[RowDim+1];
        for (int d=1; d<=RowDim; ++d) {
            Len[d] = 0;
            El[d] = NULL;
        }
        ElSorted = false;
    }

    spMatrix::~spMatrix(void)
    {
        Clear();
    }

    void spMatrix::Clear(void)
    {
        if (Len != NULL) DELARR(Len); 
        if (El != NULL) {
            for (int d=1; d<=RowDim; ++d) 
                if (El[d] != NULL) DELARR(El[d]);
            if (El != NULL) DELARR(El);
        }
    }

    void spMatrix::reshape(const size_t r, const size_t c)
    {
        int d;
   
        if (Len != NULL) DELARR(Len);
        for (d=1; d<=RowDim; ++d) 
            if (El[d] != NULL) DELARR(El[d]);
        if (El != NULL) DELARR(El);

        RowDim = r;
        ClmDim = c;
        Len = new size_t[RowDim+1];
        El = new ElType*[RowDim+1];
        for (d=1; d<=RowDim; ++d) {
            Len[d] = 0;
            El[d] = NULL;
        }
        ElSorted = false;
    }

    double spMatrix::Value(const size_t r, const size_t c) const
    {
        double Val = 0.;
   
        size_t ElCount;
        ElType *PtrEl;

        if (r+1 > 0 && r+1 <= RowDim && c+1 > 0 && c+1 <= ClmDim) {
            PtrEl = El[r+1];
            for (ElCount = Len[r+1]; ElCount > 0; ElCount--) {
                if (PtrEl->Pos == c+1)
                    Val = PtrEl->Val;
                PtrEl++;
            }
        }
        return(Val);
    }

    double& spMatrix::operator()(const size_t r, const size_t c)
    {
        size_t pos = 0;
        myBool found = FALSE;
    
        size_t ElCount;

        for (ElCount = 0; ElCount < Len[r+1] && !found; ++ElCount) {
            if (El[r+1][ElCount].Pos == c+1 || El[r+1][ElCount].Pos == 0) {
                pos = ElCount;
                El[r+1][ElCount].Pos = c+1;
                found = TRUE;
            }
        }
        if (!found) {
            ExtendRow(r, Len[r+1]+1);
            pos = Len[r+1]-1;
            El[r+1][pos].Pos = c+1;
        }
      
        return El[r+1][pos].Val;
    }

    void spMatrix::Zero(const size_t r, const size_t c)
    {
        size_t pos = 0;
        myBool found = FALSE;
    
        size_t ElCount;

        for (ElCount = 0; ElCount < Len[r+1] && !found; ++ElCount) {
            if (El[r+1][ElCount].Pos == c+1 || El[r+1][ElCount].Pos == 0) {
                pos = ElCount;
                El[r+1][ElCount].Pos = c+1;
                found = TRUE;
            }
        }
        if (found) {
            for (ElCount=0; ElCount<Len[r+1] && El[r+1][ElCount].Pos != 0; ++ElCount) ;
            El[r+1][pos] = El[r+1][ElCount-1];
            El[r+1][ElCount-1] = ZeroEl;
        }
   
    }

    void spMatrix::ExtendDim(const size_t r, const size_t c)
    {
        ClmDim = c;
        if (r > RowDim) {
            size_t* newlen = new size_t[r+1];
            ElType** newel = new ElType*[r+1];

            size_t i;
            for (i=1; i<=RowDim; ++i) {
                newlen[i] = Len[i];
                newel[i] = El[i];
            }
            for (; i<=r; ++i) {      
                newlen[i] = 0;
                newel[i] = NULL;
            }

            RowDim = r;
            DELARR(Len);
            DELARR(El);
            Len = newlen;
            El = newel;

            ElSorted = FALSE;
        }
    }
   
    void spMatrix::ExtendRow(const size_t r, const size_t e)
    {
        size_t ElCount;
        ElType *PtrEl;
        size_t ee = e == 1 ? 2 : e;

        if (r+1 > 0 && r <= RowDim) {
            PtrEl = El[r+1];

            if (ee > 0) {
                El[r+1] = new ElType[ee];
         
                for (ElCount = 0; ElCount<Len[r+1] && ElCount<e; ++ElCount) 
                    El[r+1][ElCount] = PtrEl[ElCount];

                for (; ElCount < ee; ++ElCount)
                    El[r+1][ElCount] = ZeroEl;
            } else {
                El[r+1] = NULL;
            }
      
            if (PtrEl != NULL) DELARR(PtrEl);
      
            Len[r+1] = ee;
      
        }
    }

    void spMatrix::Allocate(const size_t r, const size_t e)
    {
        size_t ElCount;
        ElType *PtrEl;

        if (r+1 > 0 && r <= RowDim) {
            Len[r+1] = e;
      
            PtrEl = El[r+1];
      
            if (PtrEl != NULL) DELARR(PtrEl);
      
            if (e > 0) {
                PtrEl = new ElType[e];
                El[r+1] = PtrEl;
         
                for (ElCount = e; ElCount > 0; ElCount--) {
                    *PtrEl = ZeroEl;
                    PtrEl++;
                }
            } else {
                El[r+1] = NULL;
            }
        }
    }

    double* spMatrix::ColSums(double* s)
    {
        int r;
        double* sum;

        if (s == NULL) {
            sum = new double[ClmDim];
            for (r=0; r<ClmDim; ++r) sum[r] = 0.;
        } else
            sum = s;

        for (r=1; r<=RowDim; ++r) 
            for (int i=0; i<Len[r]; ++i)
                if (El[r][i].Pos)
                    sum[El[r][i].Pos-1] -= El[r][i].Val;
        return sum;
    }     

    spMatrix& spMatrix::operator=(const spMatrix& A)
    {
        reshape(A.RowDim,A.ClmDim);
        int r, c;
        for (r=1; r<=RowDim; ++r) {
            Len[r] = A.Len[r];
            El[r] = new ElType[Len[r]];
            for (c=0; c<Len[r]; ++c)
                El[r][c] = A.El[r][c];
        }
        ElSorted = A.ElSorted;
        return *this;
    }


    void spMatrix::print(std::ostream& s) const
    {
        s << RowDim << 'x' << ClmDim << '\n';
        for (int i=1; i<=RowDim; ++i) {
            for (int j=0; j<Len[i]; ++j)
                if (El[i][j].Pos)
                    s << '(' << El[i][j].Pos << ',' << El[i][j].Val << ") ";
            s << '\n';
        }
    }
        
}

