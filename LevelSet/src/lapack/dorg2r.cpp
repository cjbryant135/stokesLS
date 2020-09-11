/*
 * C++ implementation of lapack routine dorg2r
 *
 * $Id: dorg2r.cpp,v 1.6 1993/04/06 20:41:38 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:39
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dorg2r.cpp,v $
 * Revision 1.6  1993/04/06 20:41:38  alv
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
 * Revision 1.2  1993/03/05  23:16:10  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:55  alv
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

RWLAPKDECL void /*FUNCTION*/ dorg2r(const long &m, const long &n, const long &k, double *a, 
                                    const long &lda, double tau[], double work[], long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, _do2, i, i_, j, j_, l, l_;

  
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
  
    //  DORG2R generates an m by n real matrix Q with orthonormal columns,
    //  which is defined as the first n columns of a product of k elementary
    //  reflectors of order m
  
    //        Q  =  H(1) H(2) . . . H(k)
  
    //  as returned by DGEQRF.
  
    //  Arguments
    //  =========
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix Q. M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix Q. M >= N >= 0.
  
    //  K       (input) INTEGER
    //          The number of elementary reflectors whose product defines the
    //          matrix Q. N >= K >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the i-th column must contain the vector which
    //          defines the elementary reflector H(i), for i = 1,2,...,k, as
    //          returned by DGEQRF in the first k columns of its array
    //          argument A.
    //          On exit, the m-by-n matrix Q.
  
    //  LDA     (input) INTEGER
    //          The first dimension of the array A. LDA >= max(1,M).
  
    //  TAU     (input) DOUBLE PRECISION array, dimension (K)
    //          TAU(i) must contain the scalar factor of the elementary
    //          reflector H(i), as returned by DGEQRF.
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
  
    //  INFO    (output) INTEGER
    //          = 0: successful exit
    //          < 0: if INFO = -i, the i-th argument has an illegal value
  
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
    else if( n < 0 || n > m ) { 
        info = -2;
    }
    else if( k < 0 || k > n ) { 
        info = -3;
    }
    else if( lda < max( 1, m ) ) { 
        info = -5;
    }
    if( info != 0 ) { 
        xerbla( "DORG2R", -info );
        return;
    }
  
    //     Quick return if possible
  
    if( n <= 0 ) 
        return;
  
    //     Initialise columns k+1:n to columns of the unit matrix
  
    for( j = k + 1, j_ = j - 1, _do0 = n; j <= _do0; j++, j_++ ) { 
        for( l = 1, l_ = l - 1, _do1 = m; l <= _do1; l++, l_++ ) { 
            A(j_,l_) = ZERO;
        }
        A(j_,j_) = ONE;
    }
  
    for( i = k, i_ = i - 1; i >= 1; i--, i_-- ) { 
    
        //        Apply H(i) to A(i:m,i:n) from the left
    
        if( i < n ) { 
            A(i_,i_) = ONE;
            dlarf( 'L'/* Left */, m - i + 1, n - i, &A(i_,i_), 1, tau[i_], 
                   &A(i_ + 1,i_), lda, work );
        }
        if( i < m ) 
            dscal( m - i, -tau[i_], &A(i_,i_ + 1), 1 );
        A(i_,i_) = ONE - tau[i_];
    
        //        Set A(1:i-1,i) to zero
    
        for( l = 1, l_ = l - 1, _do2 = i - 1; l <= _do2; l++, l_++ ) { 
            A(i_,l_) = ZERO;
        }
    }
    return;
  
    //     End of DORG2R
  
#undef  A
} // end of function 

