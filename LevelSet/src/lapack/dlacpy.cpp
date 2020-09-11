/*
 * C++ implementation of lapack routine dlacpy
 *
 * $Id: dlacpy.cpp,v 1.4 1993/04/06 20:40:53 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:34:59
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlacpy.cpp,v $
 * Revision 1.4  1993/04/06 20:40:53  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:07  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:05  alv
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

RWLAPKDECL void /*FUNCTION*/ dlacpy(const char &uplo, const long &m, const long &n, double *a, 
                                    const long &lda, double *b, const long &ldb)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
#define B(I_,J_)  (*(b+(I_)*(ldb)+(J_)))
    long _do0, _do1, _do2, _do3, _do4, _do5, i, i_, j, j_;

  
    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992
  
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  DLACPY copies all or part of a two-dimensional matrix A to another
    //  matrix B.
  
    //  Arguments
    //  =========
  
    //  UPLO    (input) CHARACTER*1
    //          Specifies the part of the matrix A to be copied to B.
    //          = 'U':      Upper triangular part
    //          = 'L':      Lower triangular part
    //          Otherwise:  All of the matrix A
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix A.  M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix A.  N >= 0.
  
    //  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
    //          The m by n matrix A.  If UPLO = 'U', only the upper triangle
    //          or trapezoid is accessed; if UPLO = 'L', only the lower
    //          triangle or trapezoid is accessed.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
    //  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
    //          On exit, B = A in the locations specified by UPLO.
  
    //  LDB     (input) INTEGER
    //          The leading dimension of the array B.  LDB >= max(1,M).
  
    //  =====================================================================
  
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    if( lsame( uplo, 'U' ) ) { 
        for( j = 1, j_ = j - 1, _do0 = n; j <= _do0; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do1 = min( j, m ); i <= _do1; i++, i_++ ) { 
                B(j_,i_) = A(j_,i_);
            }
        }
    }
    else if( lsame( uplo, 'L' ) ) { 
        for( j = 1, j_ = j - 1, _do2 = n; j <= _do2; j++, j_++ ) { 
            for( i = j, i_ = i - 1, _do3 = m; i <= _do3; i++, i_++ ) { 
                B(j_,i_) = A(j_,i_);
            }
        }
    }
    else { 
        for( j = 1, j_ = j - 1, _do4 = n; j <= _do4; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do5 = m; i <= _do5; i++, i_++ ) { 
                B(j_,i_) = A(j_,i_);
            }
        }
    }
    return;
  
    //     End of DLACPY
  
#undef  B
#undef  A
} // end of function 

