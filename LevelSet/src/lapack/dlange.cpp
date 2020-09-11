/*
 * C++ implementation of lapack routine dlange
 *
 * $Id: dlange.cpp,v 1.4 1993/04/06 20:41:03 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:35:25
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlange.cpp,v $
 * Revision 1.4  1993/04/06 20:41:03  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:24  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:17  alv
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

RWLAPKDECL double dlange(const char &norm, const long &m, const long &n, double *a, 
                         const long &lda, double work[])
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, _do2, _do3, _do4, _do5, _do6, _do7, _do8, 
        i, i_, j, j_;
    double dlange_v, scale, sum, value;

  
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
  
    //  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
    //  the  infinity norm,  or the  element of  largest absolute value  of a
    //  real matrix A.
  
    //  Description
    //  ===========
  
    //  DLANGE returns the value
  
    //     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
    //              (
    //              ( norm1(A),         NORM = '1', 'O' or 'o'
    //              (
    //              ( normI(A),         NORM = 'I' or 'i'
    //              (
    //              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
  
    //  where  norm1  denotes the  one norm of a matrix (maximum column sum),
    //  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    //  normF  denotes the  Frobenius norm of a matrix (square root of sum of
    //  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
  
    //  Arguments
    //  =========
  
    //  NORM    (input) CHARACTER*1
    //          Specifies the value to be returned in DLANGE as described
    //          above.
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix A.  M >= 0.  When M = 0,
    //          DLANGE is set to zero.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix A.  N >= 0.  When N = 0,
    //          DLANGE is set to zero.
  
    //  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
    //          The m by n matrix A.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(M,1).
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
    //          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
    //          referenced.
  
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    if( min( m, n ) == 0 ) { 
        value = ZERO;
    }
    else if( lsame( norm, 'M' ) ) { 
    
        //        Find max(abs(A(i,j))).
    
        value = ZERO;
        for( j = 1, j_ = j - 1, _do0 = n; j <= _do0; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do1 = m; i <= _do1; i++, i_++ ) { 
                value = max( value, fabs( A(j_,i_) ) );
            }
        }
    }
    else if( (lsame( norm, 'O' )) || (norm == '1') ) { 
    
        //        Find norm1(A).
    
        value = ZERO;
        for( j = 1, j_ = j - 1, _do2 = n; j <= _do2; j++, j_++ ) { 
            sum = ZERO;
            for( i = 1, i_ = i - 1, _do3 = m; i <= _do3; i++, i_++ ) { 
                sum = sum + fabs( A(j_,i_) );
            }
            value = max( value, sum );
        }
    }
    else if( lsame( norm, 'I' ) ) { 
    
        //        Find normI(A).
    
        for( i = 1, i_ = i - 1, _do4 = m; i <= _do4; i++, i_++ ) { 
            work[i_] = ZERO;
        }
        for( j = 1, j_ = j - 1, _do5 = n; j <= _do5; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do6 = m; i <= _do6; i++, i_++ ) { 
                work[i_] = work[i_] + fabs( A(j_,i_) );
            }
        }
        value = ZERO;
        for( i = 1, i_ = i - 1, _do7 = m; i <= _do7; i++, i_++ ) { 
            value = max( value, work[i_] );
        }
    }
    else if( (lsame( norm, 'F' )) || (lsame( norm, 'E' )) ) { 
    
        //        Find normF(A).
    
        scale = ZERO;
        sum = ONE;
        for( j = 1, j_ = j - 1, _do8 = n; j <= _do8; j++, j_++ ) { 
            dlassq( m, &A(j_,0), 1, scale, sum );
        }
        value = scale*sqrt( sum );
    }
  
    dlange_v = value;
    return( dlange_v );
  
    //     End of DLANGE
  
#undef  A
} // end of function 

