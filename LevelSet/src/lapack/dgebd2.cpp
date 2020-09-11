/*
 * C++ implementation of lapack routine dgebd2
 *
 * $Id: dgebd2.cpp,v 1.6 1993/04/06 20:40:24 alv Exp $
 *
 *************************************************************
 *
 * Rogue Wave Software, Inc.
 * P.O. Box 2328
 * Corvallis, OR 97339
 *
 * Copyright (C) 1993.
 * This software is subject to copyright protection under the
 * laws of the United States and other countries.
 *
 *************************************************************
 *
 * Translated from the Fortran using Cobalt Blue's FOR_C++,
 * and then massaged slightly to Rogue Wave format.
 *
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:33:51
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dgebd2.cpp,v $
 * Revision 1.6  1993/04/06 20:40:24  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.5  1993/03/19  18:41:23  alv
 * now passes chars explicitly, rather than indirection of a string, to shut up SUN warnings
 *
 * Revision 1.4  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:14:30  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:06:30  alv
 * Initial revision
 *
 */

#if 0
#include "rw/lapkdefs.h"
#include "rw/bla.h"
#include "rw/lapack.h"
#include "rw/fortran.h" /* Fortran run time library */
#else
#include "level/lapack.h"
#endif

RWLAPKDECL void /*FUNCTION*/ dgebd2(const long &m, const long &n, double *a, const long &lda, 
                                    double d[], double e[], double tauq[], double taup[], double work[], 
                                    long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ZERO = 0.0e0;
    const double ONE = 1.0e0;
// end of PARAMETER translations

    long _do0, _do1, i, i_;

  
    //  -- LAPACK routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992
  
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  DGEBD2 reduces a real general m by n matrix A to upper or lower
    //  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.
  
    //  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
  
    //  Arguments
    //  =========
  
    //  M       (input) INTEGER
    //          The number of rows in the matrix A.  M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns in the matrix A.  N >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the m by n general matrix to be reduced.
    //          On exit,
    //          if m >= n, the diagonal and the first superdiagonal are
    //            overwritten with the upper bidiagonal matrix B; the
    //            elements below the diagonal, with the array TAUQ, represent
    //            the orthogonal matrix Q as a product of elementary
    //            reflectors, and the elements above the first superdiagonal,
    //            with the array TAUP, represent the orthogonal matrix P as
    //            a product of elementary reflectors;
    //          if m < n, the diagonal and the first subdiagonal are
    //            overwritten with the lower bidiagonal matrix B; the
    //            elements below the first subdiagonal, with the array TAUQ,
    //            represent the orthogonal matrix Q as a product of
    //            elementary reflectors, and the elements above the diagonal,
    //            with the array TAUP, represent the orthogonal matrix P as
    //            a product of elementary reflectors.
    //          See Further Details.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
    //  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
    //          The diagonal elements of the bidiagonal matrix B:
    //          D(i) = A(i,i).
  
    //  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
    //          The off-diagonal elements of the bidiagonal matrix B:
    //          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
    //          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
  
    //  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
    //          The scalar factors of the elementary reflectors which
    //          represent the orthogonal matrix Q. See Further Details.
  
    //  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
    //          The scalar factors of the elementary reflectors which
    //          represent the orthogonal matrix P. See Further Details.
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (max(M,N))
  
    //  INFO    (output) INTEGER
    //          = 0: successful exit.
    //          < 0: if INFO = -i, the i-th argument had an illegal value.
  
    //  Further Details
    //  ===============
  
    //  The matrices Q and P are represented as products of elementary
    //  reflectors:
  
    //  If m >= n,
  
    //     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
  
    //  Each H(i) and G(i) has the form:
  
    //     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  
    //  where tauq and taup are real scalars, and v and u are real vectors;
    //  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
    //  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
    //  tauq is stored in TAUQ(i) and taup in TAUP(i).
  
    //  If m < n,
  
    //     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
  
    //  Each H(i) and G(i) has the form:
  
    //     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  
    //  where tauq and taup are real scalars, and v and u are real vectors;
    //  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
    //  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
    //  tauq is stored in TAUQ(i) and taup in TAUP(i).
  
    //  The contents of A on exit are illustrated by the following examples:
  
    //  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
  
    //    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
    //    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
    //    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
    //    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
    //    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
    //    (  v1  v2  v3  v4  v5 )
  
    //  where d and e denote diagonal and off-diagonal elements of B, vi
    //  denotes an element of the vector defining H(i), and ui an element of
    //  the vector defining G(i).
  
    //  =====================================================================
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test the input parameters
  
    info = 0;
    if( m < 0 ) { 
        info = -1;
    }
    else if( n < 0 ) { 
        info = -2;
    }
    else if( lda < max( 1, m ) ) { 
        info = -4;
    }
    if( info < 0 ) { 
        xerbla( "DGEBD2", -info );
        return;
    }
  
    if( m >= n ) { 
    
        //        Reduce to upper bidiagonal form
    
        for( i = 1, i_ = i - 1, _do0 = n; i <= _do0; i++, i_++ ) { 
      
            //           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
      
            dlarfg( m - i + 1, A(i_,i_), &A(i_,min( i + 1, m ) - 1), 
                    1, tauq[i_] );
            d[i_] = A(i_,i_);
            A(i_,i_) = ONE;
      
            //           Apply H(i) to A(i:m,i+1:n) from the left
      
            dlarf( 'L'/* Left */, m - i + 1, n - i, &A(i_,i_), 1, tauq[i_], 
                   &A(i_ + 1,i_), lda, work );
            A(i_,i_) = d[i_];
      
            if( i < n ) { 
        
                //              Generate elementary reflector G(i) to annihilate
                //              A(i,i+2:n)
        
                dlarfg( n - i, A(i_ + 1,i_), &A(min( i + 2, n ) - 1,i_), 
                        lda, taup[i_] );
                e[i_] = A(i_ + 1,i_);
                A(i_ + 1,i_) = ONE;
        
                //              Apply G(i) to A(i+1:m,i+1:n) from the right
        
                dlarf( 'R'/* Right */, m - i, n - i, &A(i_ + 1,i_), 
                       lda, taup[i_], &A(i_ + 1,i_ + 1), lda, work );
                A(i_ + 1,i_) = e[i_];
            }
            else { 
                taup[i_] = ZERO;
            }
        }
    }
    else { 
    
        //        Reduce to lower bidiagonal form
    
        for( i = 1, i_ = i - 1, _do1 = m; i <= _do1; i++, i_++ ) { 
      
            //           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
      
            dlarfg( n - i + 1, A(i_,i_), &A(min( i + 1, n ) - 1,i_), 
                    lda, taup[i_] );
            d[i_] = A(i_,i_);
            A(i_,i_) = ONE;
      
            //           Apply G(i) to A(i+1:m,i:n) from the right
      
            dlarf( 'R'/* Right */, m - i, n - i + 1, &A(i_,i_), lda, 
                   taup[i_], &A(i_,min( i + 1, m ) - 1), lda, work );
            A(i_,i_) = d[i_];
      
            if( i < m ) { 
        
                //              Generate elementary reflector H(i) to annihilate
                //              A(i+2:m,i)
        
                dlarfg( m - i, A(i_,i_ + 1), &A(i_,min( i + 2, m ) - 1), 
                        1, tauq[i_] );
                e[i_] = A(i_,i_ + 1);
                A(i_,i_ + 1) = ONE;
        
                //              Apply H(i) to A(i+1:m,i+1:n) from the left
        
                dlarf( 'L'/* Left */, m - i, n - i, &A(i_,i_ + 1), 
                       1, tauq[i_], &A(i_ + 1,i_ + 1), lda, work );
                A(i_,i_ + 1) = e[i_];
            }
            else { 
                tauq[i_] = ZERO;
            }
        }
    }
    return;
  
    //     End of DGEBD2
  
#undef  A
} // end of function 

