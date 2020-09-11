#include <iostream>
#include "spqmatrix.h"

namespace blockmatrix {
        
    spQMatrix::spQMatrix(size_t d)
    /* constructor of the type QMatrix */
    {
        size_t RoC;

        Dim = d;
        Len = new size_t[Dim+1];
        El = new ElType*[Dim+1];
        DiagEl = new ElType*[Dim+1];
        InvDiagEl = new double[Dim+1];
//   ILU = (QMatrix *)malloc(sizeof(QMatrix));
        for (RoC = 1; RoC <= Dim; RoC++) {
            Len[RoC] = 0;
            El[RoC] = NULL;
            DiagEl[RoC] = NULL;
            InvDiagEl[RoC] = 0.0;
        }
        ElSorted = FALSE;
        DiagElAlloc = FALSE;
        ZeroInDiag = TRUE;
    }

    spQMatrix::~spQMatrix(void)
    {
        Clear();
    }

    void spQMatrix::Clear(void)
    {
        size_t RoC;

        if (Len != NULL && El != NULL) {
            for (RoC = 1; RoC <= Dim; RoC++) {
                if (Len[RoC] > 0) {
                    if (El[RoC] != NULL)
                        delete[] El[RoC];
                }
            }
        }
        if (Len != NULL) {
            delete[] Len;
            Len = NULL;
        }
        if (El != NULL) {
            delete[] El;
            El = NULL;
        }
        if (DiagEl != NULL) {
            delete[] DiagEl;
            DiagEl = NULL;
        }
        if (InvDiagEl != NULL) {
            delete[] InvDiagEl;
            InvDiagEl = NULL;
        }
    }

    void spQMatrix::reshape(const size_t r, const size_t c)
    {
        size_t RoC;

        if (Len != NULL && El != NULL) {
            for (RoC = 1; RoC <= Dim; RoC++) {
                if (Len[RoC] > 0) {
                    if (El[RoC] != NULL)
                        delete[] El[RoC];
                }
            }
        }
        if (Len != NULL) {
            delete[] Len;
            Len = NULL;
        }
        if (El != NULL) {
            delete[] El;
            El = NULL;
        }
        if (DiagEl != NULL) {
            delete[] DiagEl;
            DiagEl = NULL;
        }
        if (InvDiagEl != NULL) {
            delete[] InvDiagEl;
            InvDiagEl = NULL;
        }

        Dim = r;
        Len = new size_t[Dim+1];
        El = new ElType*[Dim+1];
        DiagEl = new ElType*[Dim+1];
        InvDiagEl = new double[Dim+1];
//   ILU = (QMatrix *)malloc(sizeof(QMatrix));
        for (RoC = 1; RoC <= Dim; RoC++) {
            Len[RoC] = 0;
            El[RoC] = NULL;
            DiagEl[RoC] = NULL;
            InvDiagEl[RoC] = 0.0;
        }
        ElSorted = FALSE;
        DiagElAlloc = FALSE;
        ZeroInDiag = TRUE;
    }

    double spQMatrix::Value(const size_t r, const size_t c) const
    {
        double Val;
   
        size_t ElCount;
        ElType *PtrEl;

        if (r+1 > 0 && r+1 <= Dim && c+1 > 0 && c+1 <= Dim) {
            Val = 0.0;
            PtrEl = El[r+1];
            for (ElCount = Len[r+1]; ElCount > 0; ElCount--) {
                if (PtrEl->Pos == c+1)
                    Val = PtrEl->Val;
                PtrEl++;
            }

        }
        return(Val);
    }

    double& spQMatrix::operator()(const size_t r, const size_t c)
    {
        size_t pos = 0;
        myBool found = FALSE;
    
        size_t ElCount;

        if (r+1 > Dim) 
            ExtendDim(r+1);

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

    void spQMatrix::Zero(const size_t r, const size_t c)
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

    void spQMatrix::Allocate(const size_t r, const size_t e)
    {
        size_t ElCount;
        ElType *PtrEl;

        if (r+1 > 0 && r <= Dim) {
            Len[r+1] = e;
      
            PtrEl = El[r+1];
      
            if (PtrEl != NULL) {
                delete[] PtrEl;
                PtrEl = NULL;
            }
      
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

    void spQMatrix::ExtendDim(const size_t r)
    {
        if (r > Dim) {
            size_t* newlen = new size_t[r+1];
            ElType** newel = new ElType*[r+1];
            ElType** newdiag = new ElType*[r+1];
            double* newidiag = new double[r+1];

            size_t i;
            for (i=1; i<=Dim; ++i) {
                newlen[i] = Len[i];
                newel[i] = El[i];
                newdiag[i] = DiagEl[i];
                newidiag[i] = InvDiagEl[i];
            }
            for (; i<=r; ++i) {      
                newlen[i] = 0;
                newel[i] = NULL;
                newdiag[i] = NULL;
                newidiag[i] = 0.;
            }

            Dim = r;
            delete[] Len;
            delete[] El;
            delete[] DiagEl;
            delete[] InvDiagEl;
            Len = newlen;
            El = newel;
            DiagEl = newdiag;
            InvDiagEl = newidiag;

            ElSorted = FALSE;
            DiagElAlloc = FALSE;
            ZeroInDiag = TRUE;
        }
    }
   
    void spQMatrix::ExtendRow(const size_t r, const size_t e)
    {
        size_t ElCount;
        ElType *PtrEl;
        size_t ee = e == 1 ? 2 : e;

        if (r+1 > 0 && r <= Dim) {
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
      
            if (PtrEl != NULL) {
                delete[] PtrEl;
                PtrEl = NULL;
            }
      
            Len[r+1] = ee;
      
        }
    }

    static int ElCompar(const void *El1, const void *El2)
    /* compares positions of two matrix elements */
    {
        int Compar;

        Compar = 0;
        if (((ElType *)El1)->Pos < ((ElType *)El2)->Pos)
            Compar = -1;
        if (((ElType *)El1)->Pos > ((ElType *)El2)->Pos)
            Compar = +1;

        return(Compar);
    }

    void spQMatrix::SortEl(void)
    {
        size_t RoC;
        myBool UpperOnly;
   
        UpperOnly = TRUE;
        for (RoC = 1; RoC <= Dim; RoC++) 
            qsort((void *)El[RoC], Len[RoC], sizeof(ElType), ElCompar);
      
        ElSorted = TRUE;
        DiagElAlloc = FALSE;
        ZeroInDiag = TRUE;  
    }

    void spQMatrix::AllocInvDiagEl(void)
    /* allocate pointers and compute inverse for diagonal elements of the matrix */
    {
        size_t RoC, ElCount;
        myBool Found;
        ElType *PtrEl;
   
        if (!DiagElAlloc) {
            ZeroInDiag = FALSE;
            for (RoC = 1; RoC <= Dim; RoC++) {
                Found = FALSE;
                PtrEl = El[RoC] + Len[RoC] - 1;
                for (ElCount = Len[RoC]; ElCount > 0; ElCount--) {
                    if (PtrEl->Pos == RoC) {
                        Found = TRUE;
                        DiagEl[RoC] = PtrEl;
                    }
                    PtrEl--;
                }
                if (!Found) {
                    ZeroInDiag = TRUE;
                    DiagEl[RoC] = (Element*)(&ZeroEl);
                }
            }
            DiagElAlloc = TRUE;
        
            if (!ZeroInDiag) {
                for (RoC = 1; RoC <= Dim; RoC++)
                    InvDiagEl[RoC] = 1.0 / DiagEl[RoC]->Val;
            }
        }
    }

    void spQMatrix::SwapRows(const size_t r1, const size_t r2)
    {
        ElType *TempRow = El[r1+1];
        El[r1+1] = El[r2+1];
        El[r2+1] = TempRow;
        size_t temp = Len[r1+1];
        Len[r1+1] = Len[r2+1];
        Len[r2+1] = temp;
    }

    double* spQMatrix::ColSums(double* s)
    {
        int r;
        double* sum;

        if (s == NULL) {
            sum = new double[Dim];
            for (r=0; r<Dim; ++r) sum[r] = 0.;
        } else
            sum = s;
   
        for (r=1; r<=Dim; ++r) 
            for (int i=0; i<Len[r]; ++i)
                if (El[r][i].Pos)
                    sum[El[r][i].Pos-1] += El[r][i].Val;
        return sum;
    }     

    spQMatrix& spQMatrix::operator=(const spQMatrix& A)
    {
        reshape(A.Dim);
        int r, c;
        for (r=1; r<=Dim; ++r) {
            Len[r] = A.Len[r];
            El[r] = new ElType[Len[r]];
            for (c=0; c<Len[r]; ++c)
                El[r][c] = A.El[r][c];
        }
        ElSorted = A.ElSorted;
        DiagElAlloc = FALSE;
        ZeroInDiag = A.ZeroInDiag;
        return *this;
    }

    void spQMatrix::print(std::ostream& s) const
    {
        s << Dim << 'x' << Dim << '\n';
        for (int i=1; i<=Dim; ++i) {
            for (int j=0; j<Len[i]; ++j)
                if (El[i][j].Pos)
                    s << '(' << El[i][j].Pos << ',' << El[i][j].Val << ") ";
            s << '\n';
        }
    }

    void spQMatrix::matlab_print(const char* n, std::ostream& s) const
    {
        s << n << " = zeros(" << Dim << ',' << Dim << ");\n";
        for (int i=1; i<=Dim; ++i) {
            for (int j=0; j<Len[i]; ++j)
                if (El[i][j].Pos)
                    s << n << "(" << i << ',' << El[i][j].Pos << ") = " << El[i][j].Val << ";\n";
        }
        s.flush();
    }
        
}
