/*
 * C++ implementation of lapack routine dgelq2
 *
 * $Id: dgelq2.cpp,v 1.6 1993/04/06 20:40:30 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:34:07
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dgelq2.cpp,v $
 * Revision 1.6  1993/04/06 20:40:30  alv
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
 * Revision 1.2  1993/03/05  23:14:38  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:06:37  alv
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

RWLAPKDECL void /*FUNCTION*/ dgelq2(const long &m, const long &n, double *a, const long &lda, 
                                    double tau[], double work[], long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
// end of PARAMETER translations

    long _do0, i, i_, k;
    double aii;

  
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
  
    //  DGELQ2 computes an LQ factorization of a real m by n matrix A:
    //  A = L * Q.
  
    //  Arguments
    //  =========
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix A.  M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix A.  N >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the m by n matrix A.
    //          On exit, the elements on and below the diagonal of the array
    //          contain the m by min(m,n) lower trapezoidal matrix L (L is
    //          lower triangular if m <= n); the elements above the diagonal,
    //          with the array TAU, represent the orthogonal matrix Q as a
    //          product of elementary reflectors (see Further Details).
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
    //  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
    //          The scalar factors of the elementary reflectors (see Further
    //          Details).
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
  
    //  INFO    (output) INTEGER
    //          = 0: successful exit
    //          < 0: if INFO = -i, the i-th argument had an illegal value
  
    //  Further Details
    //  ===============
  
    //  The matrix Q is represented as a product of elementary reflectors
  
    //     Q = H(k) . . . H(2) H(1), where k = min(m,n).
  
    //  Each H(i) has the form
  
    //     H(i) = I - tau * v * v'
  
    //  where tau is a real scalar, and v is a real vector with
    //  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
    //  and tau in TAU(i).
  
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
  
    //     Test the input arguments
  
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
    if( info != 0 ) { 
        xerbla( "DGELQ2", -info );
        return;
    }
  
    k = min( m, n );
  
    for( i = 1, i_ = i - 1, _do0 = k; i <= _do0; i++, i_++ ) { 
    
        //        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
    
        dlarfg( n - i + 1, A(i_,i_), &A(min( i + 1, n ) - 1,i_), lda, 
                tau[i_] );
        if( i < m ) { 
      
            //           Apply H(i) to A(i+1:m,i:n) from the right
      
            aii = A(i_,i_);
            A(i_,i_) = ONE;
            dlarf( 'R'/* Right */, m - i, n - i + 1, &A(i_,i_), lda, 
                   tau[i_], &A(i_,i_ + 1), lda, work );
            A(i_,i_) = aii;
        }
    }
    return;
  
    //     End of DGELQ2
  
#undef  A
} // end of function 

