/*************************************************
    spglobals.cpp

    $Header: spglobals.cpp,v 1.1 99/09/20 11:35:51 chopp Exp $

    $Log:       spglobals.cpp,v $
Revision 1.1  99/09/20  11:35:51  11:35:51  chopp (David Chopp)
Initial revision

Revision 1.4  99/04/05  13:32:20  13:32:20  chopp (David Chopp)
*** none ***

Revision 1.3  99/03/09  10:12:44  10:12:44  chopp (David Chopp)
*** none ***

Revision 1.1  99/02/26  14:02:39  14:02:39  chopp (David Chopp)
Initial revision

Revision 1.1  99/02/26  14:01:52  14:01:52  chopp (David Chopp)
Initial revision

*************************************************/

#include "spmatrix.h"
#include "spqmatrix.h"
#include "spvector.h"
#include <iostream>
#include <stdlib.h>
#ifdef PLOT_RHO
#include "plotwindow3d.h"
extern double* com[2];
#endif

namespace blockmatrix {
        
    spVector operator*(const spQMatrix& Q, const spVector& V)
    {
        spVector VRes(V.length());

        double Sum;
        size_t Dim, Row, Len, ElCount;
        ElType *PtrEl;

        Dim = V.length();

        /* multiplication of the lower, diagonal and upper part
           of the matrix Q by the vector V */
        for (Row = 1; Row <= Dim; Row++) {
            Len = Q.Len[Row];
            PtrEl = Q.El[Row];
            Sum = 0.0;
            for (ElCount = Len; ElCount > 0; ElCount--) {
                Sum += PtrEl->Val * V.Cmp[PtrEl->Pos];
                PtrEl++;
            }
            VRes.Cmp[Row] = Sum;
        }
            
        return(VRes);
    }

    spVector operator*(const spMatrix& M, const spVector& V)
    {
        spVector VRes(M.rows());

        double Sum;
        size_t RowDim, ClmDim, Row, Len, ElCount;
        ElType *PtrEl;
    
        RowDim = M.rows();
        ClmDim = M.cols();

        /* multiplication of matrix elements by vector components */
        for (Row = 1; Row <= RowDim; Row++) {
            Len = M.Len[Row];
            PtrEl = M.El[Row];
            Sum = 0.0;
            for (ElCount = Len; ElCount > 0; ElCount--) {
                Sum += PtrEl->Val * V.Cmp[PtrEl->Pos];
                PtrEl++;
            }
            VRes.Cmp[Row] = Sum;
        }
            
        return(VRes);
    }

    double operator*(const spVector& V, const spVector& W)
    {
        double dot = 0.;
        for (int d=1; d<=V.length(); ++d)
            dot += V.Cmp[d]*W.Cmp[d];
        return dot;
    }

#if defined(__HP_aCC)
    int finite(const double x) 
    {
        return fpclassify(x) < 3;
    }
#endif

    spVector Unknown_Solve(spQMatrix& A, const spVector& b, const int MaxIter,
                           const double Omega, const double Eps)
    {
        int Iter;
        double bNorm;
        spVector r(b.length());
        spVector p(b.length());
        spVector x(b.length());
        A.SortEl();
        A.AllocInvDiagEl();
        bNorm = b.l2Norm();

        Iter = 0;
        /* x(0) = b */
        x = b;
        /* r(0) = A * b */
        p = A*b;
        r = p;
        r -= x;
        r += b;
        double rNorm = r.l2Norm();
   
        while (rNorm >= Eps * bNorm && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
               && (Iter < MaxIter || MaxIter == -1)) {
            Iter++;
            /* x(i+1) = x(i) + r(i) */
            x += r;
            if (Iter < MaxIter || MaxIter == -1) {
                /* r = b - (I - A) * x(i) */
                r = A*x;
                r -= x;
                r += b;
                rNorm = r.l2Norm();
                /* p(i+1) = A^(i+1)*b */
                if (!finite(rNorm)) {
                    dump_mat(A, b, 0);
                }
            }
        }

        std::cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        return(x);
    }

    spVector Jacobi_Solve(spQMatrix& A, const spVector& b, const int MaxIter,
                          const double Omega, const double Eps)
    {
        int Iter;
        double bNorm;
        spVector r(b.length());
        static spVector x(b.length());
        A.SortEl();
        A.AllocInvDiagEl();
        bNorm = b.l2Norm();

        Iter = 0;
        /* r = b - A * x(i) */
        if (x.length() != b.length()) 
            x.resize(b.length());
        r = b;
        r -= A*x;
        double rNorm = r.l2Norm();
   
        while (rNorm >= Eps * bNorm && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
               && Iter < MaxIter) {
            Iter++;
            /* x(i+1) = x(i) + Omega * D^(-1) * r */
            r *= A.InvDiagEl;
            r *= Omega;
            x += r;
            if (Iter < MaxIter) {
                r = b;
                r -= A*x;
            }
            rNorm = r.l2Norm();
        }

        std::cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        return(x);
    }

    spVector CG_Solve(spQMatrix& A, const spVector& b, const int MaxIter,
                      const double Omega, const double Eps)
    {
        int Iter;
        double Alpha, Beta, Rho, RhoOld = 0.0;
        double bNorm;
        spVector r(b.length());
        spVector p(b.length());
        spVector q(b.length());
        spVector z(b.length());
        spVector x(b.length());

        x = 0.;
        bNorm = b.l2Norm();
        
        Iter = 0;
        /* r = b - A * x(i) */
        if (!IsZero(x.l1Norm() / b.length())) {
            r = b;
            r -= A*x;
        } else {
            r = b;
        }
        /* plain CG (z = r) */
        double rNorm = r.l2Norm();
//   while (rNorm >= Eps * bNorm && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
        while (rNorm >= Eps && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
               && (Iter < MaxIter || MaxIter == -1)) {
            Iter++;
            Rho = r.l2Norm();
            if (Iter == 1) {
                p = r;
            } else {
                Beta = Rho * Rho / RhoOld / RhoOld;
                p *= Beta/Alpha;
                p += r;
            }
            q = A*p;
            Alpha = Rho * Rho / (p * q);
            p *= Alpha;
            x += p;
            q *= Alpha;
            r -= q;
            RhoOld = Rho;
            rNorm = r.l2Norm();
            if (!finite(rNorm)) {
                dump_mat(A, b, 0);
            }
        }

        std::cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        if (rNorm > 1.0e+5) dump_mat(A,b,0);
        return(x);
    }


    spVector SOR_Solve(spQMatrix& A, const spVector& b, const int MaxIter,
                       const double Omega, const double Eps)
    {
        int Iter;
        double bNorm;
        spVector r(b.length());
        static spVector x(b.length());
        A.SortEl();
        bNorm = b.l2Norm();

        Iter = 0;
        /* r = b - A * x(i) */
        if (x.length() != b.length()) 
            x.resize(b.length());
        x = 0.;
        r = b;
        r -= A*x;
        double rNorm = r.l2Norm();
#ifdef PLOT_RHO
        PlotWindow3D disp(true, false);
        disp.SetRange(0., 500., 0., 500., 0., 5.);
        disp.SetView(-5000,-10000,2500);
#endif
   
//   while (rNorm >= Eps * bNorm && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
        while (rNorm >= Eps && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
               && (Iter < MaxIter || MaxIter == -1)) {
            Iter++;
            for (size_t Row=1; Row <= A.size(); ++Row) 
                for (size_t e=0; e<A.Len[Row]; ++e) {
                    if (A.El[Row][e].Pos < Row) 
                        r.Cmp[Row] -= A.El[Row][e].Val*r.Cmp[A.El[Row][e].Pos];
                    else if (A.El[Row][e].Pos == Row)
                        r.Cmp[Row] /= A.El[Row][e].Val;
                }
            r *= Omega;
            x += r;
#ifdef PLOT_RHO
            if (Iter%100 == 0) {
                disp.Clear();
                for (size_t Cell=1; Cell <= x.length(); ++Cell) 
                    disp.Dot(com[0][Cell-1], com[1][Cell-1], 250.*x.Cmp[Cell]);
            }
#endif
            if (Iter < MaxIter || MaxIter == -1) {
                r = b;
                r -= A*x;
            }
            rNorm = r.l2Norm();
            if (!finite(rNorm)) {
                dump_mat(A, b, 1);
            }
//      cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        }

#ifdef PLOT_RHO
        disp.Clear();
        for (size_t Cell=1; Cell <= x.length(); ++Cell)
            disp.Dot(com[0][Cell-1], com[1][Cell-1], 250.*x.Cmp[Cell]);
#endif

        std::cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        return(x);
    }

    spVector Direct_Solve(spQMatrix& B, const spVector& b, const int MaxIter,
                          const double Omega, const double Eps)
    {
        spQMatrix A(2);
        A = B;
        size_t r, rr;
        size_t ElCount;
        spVector x(b.length());
        A.SortEl();

        x = b;
        // Forward elimination
        for (r=1; r<A.Dim; ++r) {
            if (A.Value(r-1,r-1) == 0.) {
                // zero pivot, find a non-zero pivot
                for (rr=r+1; rr<=A.Dim && A.Value(rr-1,r-1)==0.; ++rr) ;
                A.SwapRows(r-1,rr-1);
                double temp = x(r-1);
                x(r-1) = x(rr-1);
                x(rr-1) = temp;
            }
            for (rr=r+1; rr<=A.Dim; ++rr) {
                if (A.Value(rr-1,r-1) != 0.) {
                    double f = -A.Value(rr-1,r-1)/A.Value(r-1,r-1);
                    for (ElCount=0; ElCount < A.Len[r]; ++ElCount) {
                        if (A.El[r][ElCount].Pos != 0 && A.El[r][ElCount].Pos != r) 
                            A(rr-1,A.El[r][ElCount].Pos-1) += f*A.El[r][ElCount].Val;
                    }
                    A.Zero(rr-1,r-1);
                    x(rr-1) += f*x(r-1);
                }
            }
        }
   
        // Backward substitution
        for (r=A.Dim; r>1; --r) {
            x(r-1) /= A.Value(r-1,r-1);
            for (rr=r-1; rr>=1; --rr) 
                if (A.Value(rr-1,r-1) != 0.) 
                    x(rr-1) -= A.Value(rr-1,r-1)*x(r-1);
        }
                        
        return(x);
    }

    void dump_mat(spQMatrix& A, const spVector& b, const int fix)
    {
        std::cout << "A = \n";
        A.print(std::cout);
        std::cout << "\nb = \n";
        b.print(std::cout);
        spVector x(b.length());
   
        if (fix) 
            for (int i=1; i<=A.Dim; ++i) 
                for (int j=0; j<A.Len[i]; ++j) 
                    if (A.El[i][j].Pos) {
                        if (A.El[i][j].Pos == i) 
                            A.El[i][j].Val = 1-A.El[i][j].Val;
                        else
                            A.El[i][j].Val *= -1;
                    }

        x = b;
        for (int i=0; i<1000; ++i) {
            x = A*x;
            x /= x.l2Norm();
        }
        std::cout << "\nx = \n";
        x.print(std::cout);
        int l;
        double mx = 0;
        for (int k=1; k<=x.length(); ++k) 
            if (fabs(x.Cmp[k]) > mx) {
                mx = fabs(x.Cmp[k]);
                l = k;
            }
        std::cout << "\n largest index = " << l << '\n';
        exit(1);
    }

#if 0
    spVector Laplacian(spQMatrix& A, const spVector& b, const int MaxIter,
                       const double Dt, const double Eps)
    {
        int Iter;
        double bNorm;
        spVector r(b.length());
        spVector p(b.length());
        spVector x(b.length());
        A.SortEl();
        A.AllocInvDiagEl();
        bNorm = b.l2Norm();

        Iter = 0;
        /* x(0) = b */
        x = b;
        /* r(0) = A * b */
        r = A*x;
        r *= Dt;
        r += b;
        r -= x;
        double rNorm = r.l2Norm();
   
        while (rNorm >= Eps * bNorm && (!IsZero(bNorm) || !IsOne(1.0 + rNorm))
               && (Iter < MaxIter || MaxIter == -1)) {
            Iter++;
            /* x(i+1) = x(i) + r(i) */
            x += r;
            if (Iter < MaxIter || MaxIter == -1) {
                /* r = b - (I - A) * x(i) */
                r = A*x;
                r += b;
                r -= x;
                rNorm = r.l2Norm();
                /* p(i+1) = A^(i+1)*b */
                if (!finite(rNorm)) {
                    dump_mat(A, b, 0);
                }
            }
        }

        cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        return(x);
    }
#endif

}
