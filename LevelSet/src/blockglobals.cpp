#include "blockqmatrix.h"
#include "blockvector.h"
#include <iostream>
#include <stdlib.h>
#ifdef PLOT_RHO
#include "plotwindow3d.h"
#endif

namespace blockmatrix {
        
    extern double* com[2];

    BlockVector operator*(const BlockQMatrix& Q, const BlockVector& V)
    {
        BlockVector U(V.Dim);
        
        for (size_t i=0; i<V.Dim; ++i) {
            spVector* subu = new spVector(V.vec[i]->Dim);
            U.SubVector(i, *subu);
        }
        U = 0.;
        
        for (size_t i=0; i<Q.Dim; ++i) 
            for (size_t j=0; j<Q.Dim; ++j) {
                if (i==j)
                    *U.vec[i] += (*Q.mat[i][j].spq) * (*V.vec[j]);
                else
                    *U.vec[i] += (*Q.mat[i][j].sp) * (*V.vec[j]);
            }
        
        return U;
    }

    double operator*(const BlockVector& V, const BlockVector& W)
    {
        double dot = 0.;
        for (int d=0; d<V.Dim; ++d)
            dot += (*V.vec[d])*(*W.vec[d]);
        return dot;
    }

#if defined(__HP_aCC)
    int finite(const double x) 
    {
        return fpclassify(x) < 3;
    }
#endif

    BlockVector Unknown_Solve(BlockQMatrix& A, const BlockVector& b, const int MaxIter,
                              const double Omega, const double Eps)
    {
        int Iter;
        double bNorm;
        BlockVector x = b;
        BlockVector p = A*b;
        BlockVector r = p;
//   A.SortEl();
//   A.AllocInvDiagEl();
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
//      if (!finite(rNorm)) {
//            dump_mat(A, b, 0);
//         }
            }
        }

        std::cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        return(x);
    }

    BlockVector Jacobi_Solve(BlockQMatrix& A, const BlockVector& b, const int MaxIter,
                             const double Omega, const double Eps)
    {
        int Iter;
        double bNorm;
        BlockVector r = b;
        BlockVector x = b;
        A.SortEl();
        A.AllocInvDiagEl();
        bNorm = b.l2Norm();

        Iter = 0;
        /* r = b - A * x(i) */
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

    BlockVector CG_Solve(BlockQMatrix& A, const BlockVector& b, const int MaxIter,
                         const double Omega, const double Eps)
    {
        int Iter;
        double Alpha, Beta, Rho, RhoOld = 0.0;
        double bNorm;
        BlockVector r = b;
        BlockVector p = b;
        BlockVector q = b;
        BlockVector z = b;
        BlockVector x = b;

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
//      if (!finite(rNorm)) {
//         dump_mat(A, b, 0);
//      }
        }

        std::cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
//   if (rNorm > 1.0e+5) dump_mat(A,b,0);
        return(x);
    }

#if 0
    BlockVector SOR_Solve(BlockQMatrix& A, const BlockVector& b, const int MaxIter,
                          const double Omega, const double Eps)
    {
        int Iter;
        double bNorm;
        spVector r=b;
        static spVector x=b;
        A.SortEl();
        bNorm = b.l2Norm();

        Iter = 0;
        /* r = b - A * x(i) */
        x = 0.;
        r -= A*x;
        double rNorm = r.l2Norm();
#ifdef PLOT_RHO
        PlotWindow3D disp(true, false);
        disp.SetRange(0., 500., 0., 500., 0., 5.);
        disp.SetView(-5000,-10000,2500);
#endif
   
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

        cout << "Iter = " << Iter << " and rNorm = " << rNorm << '\n';
        return(x);
    }
#endif

#if 0
    spVector Direct_Solve(spQMatrix& B, const spVector& b, const int MaxIter,
                          const double Omega, const double Eps)
    {
        spQMatrix A(2);
        A = B;
        size_t r, rr, c;
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
#endif

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
