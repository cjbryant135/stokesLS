/*
 * C++ implementation of lapack routine dorgbr
 *
 * $Id: dorgbr.cpp,v 1.5 1993/04/06 20:41:38 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:40
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dorgbr.cpp,v $
 * Revision 1.5  1993/04/06 20:41:38  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:16:11  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:56  alv
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

RWLAPKDECL void /*FUNCTION*/ dorgbr(const char &vect, const long &m, const long &n, const long &k, 
                                    double *a, const long &lda, double tau[], double work[], const long &lwork, 
                                    long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ZERO = 0.0e0;
    const double ONE = 1.0e0;
// end of PARAMETER translations

    int wantq;
    long _do0, _do1, _do2, _do3, i, i_, iinfo, j, j_;

  
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
  
    //  DORGBR generates one of the matrices Q or P' determined by DGEBRD
    //  when reducing a real matrix A to bidiagonal form: A = Q * B * P'.
    //  Q and P' are defined as products of elementary reflectors H(i) or
    //  G(i) respectively..
  
    //  If VECT = 'Q', A is assumed to have been an m-by-k matrix, and Q
    //  is of order m:
    //  if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
    //  columns of Q, where m >= n >= k;
    //  if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
    //  m-by-m matrix.
  
    //  If VECT = 'P', A is assumed to have been a k-by-n matrix, and P'
    //  is of order n:
    //  if k < n, P' = G(k) . . . G(2) G(1) and DORGBR returns the first m
    //  rows of P', where n >= m >= k;
    //  if k >= n, P' = G(n-1) . . . G(2) G(1) and DORGBR returns P' as an
    //  n-by-n matrix.
  
    //  Arguments
    //  =========
  
    //  VECT    (input) CHARACTER*1
    //          Specifies whether the matrix Q or the matrix P' is required,
    //          as defined in the transformation applied by DGEBRD:
    //          = 'Q': generate Q;
    //          = 'P': generate P'.
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix Q or P' to be returned.
    //          M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix Q or P' to be returned.
    //          N >= 0.
    //          If VECT = 'Q', M >= N >= min(M,K);
    //          if VECT = 'P', N >= M >= min(N,K).
  
    //  K       (input) INTEGER
    //          If VECT = 'Q', the number of columns in the original m-by-k
    //          matrix reduced by DGEBRD.
    //          If VECT = 'P', the number of rows in the original k-by-n
    //          matrix reduced by DGEBRD.
    //          K >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the vectors which define the elementary reflectors,
    //          as returned by DGEBRD.
    //          On exit, the m-by-n matrix Q or P'.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A. LDA >= max(1,M).
  
    //  TAU     (input) DOUBLE PRECISION array, dimension
    //                                (min(M,K)) if VECT = 'Q'
    //                                (min(N,K)) if VECT = 'P'
    //          TAU(i) must contain the scalar factor of the elementary
    //          reflector H(i) or G(i), which determines Q or P', as returned
    //          by DGEBRD in its array argument TAUQ or TAUP.
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
    //          On exit, if INFO = 0, WORK(1) returns the minimum value of
    //          LWORK required to use the optimal blocksize.
  
    //  LWORK   (input) INTEGER
    //          The dimension of the array WORK. LWORK >= min(M,N).
    //          For optimum performance LWORK should be at least min(M,N)*NB,
    //          where NB is the optimal blocksize.
  
    //  INFO    (output) INTEGER
    //          = 0: successful exit
    //          < 0: if INFO = -i, the i-th argument had an illegal value
  
    //  =====================================================================
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test the input arguments
  
    info = 0;
    wantq = lsame( vect, 'Q' );
    if( !wantq && !lsame( vect, 'P' ) ) { 
        info = -1;
    }
    else if( m < 0 ) { 
        info = -2;
    }
    else if( (n < 0 || (wantq && (n > m || n < min( m, k )))) || (!wantq
                                                                  && (m > n || m < min( n, k ))) ) { 
        info = -3;
    }
    else if( k < 0 ) { 
        info = -4;
    }
    else if( lda < max( 1, m ) ) { 
        info = -6;
    }
    else if( lwork < max( 1, min( m, n ) ) ) { 
        info = -9;
    }
    if( info != 0 ) { 
        xerbla( "DORGBR", -info );
        return;
    }
  
    //     Quick return if possible
  
    if( m == 0 || n == 0 ) { 
        work[0] = 1;
        return;
    }
  
    if( wantq ) { 
    
        //        Form Q, determined by a call to DGEBRD to reduce an m-by-k
        //        matrix
    
        if( m >= k ) { 
      
            //           If m >= k, assume m >= n >= k
      
            dorgqr( m, n, k, a, lda, tau, work, lwork, iinfo );
      
        }
        else { 
      
            //           If m < k, assume m = n
      
            //           Shift the vectors which define the elementary reflectors one
            //           column to the right, and set the first row and column of Q
            //           to those of the unit matrix
      
            for( j = m, j_ = j - 1; j >= 2; j--, j_-- ) { 
                A(j_,0) = ZERO;
                for( i = j + 1, i_ = i - 1, _do0 = m; i <= _do0; i++, i_++ ) { 
                    A(j_,i_) = A(j_ - 1,i_);
                }
            }
            A(0,0) = ONE;
            for( i = 2, i_ = i - 1, _do1 = m; i <= _do1; i++, i_++ ) { 
                A(0,i_) = ZERO;
            }
            if( m > 1 ) { 
        
                //              Form Q(2:m,2:m)
        
                dorgqr( m - 1, m - 1, m - 1, &A(1,1), lda, tau, work, 
                        lwork, iinfo );
            }
        }
    }
    else { 
    
        //        Form P', determined by a call to DGEBRD to reduce a k-by-n
        //        matrix
    
        if( k < n ) { 
      
            //           If k < n, assume k <= m <= n
      
            dorglq( m, n, k, a, lda, tau, work, lwork, iinfo );
      
        }
        else { 
      
            //           If k >= n, assume m = n
      
            //           Shift the vectors which define the elementary reflectors one
            //           row downward, and set the first row and column of P' to
            //           those of the unit matrix
      
            A(0,0) = ONE;
            for( i = 2, i_ = i - 1, _do2 = n; i <= _do2; i++, i_++ ) { 
                A(0,i_) = ZERO;
            }
            for( j = 2, j_ = j - 1, _do3 = n; j <= _do3; j++, j_++ ) { 
                for( i = j - 1, i_ = i - 1; i >= 2; i--, i_-- ) { 
                    A(j_,i_) = A(j_,i_ - 1);
                }
                A(j_,0) = ZERO;
            }
            if( n > 1 ) { 
        
                //              Form P'(2:n,2:n)
        
                dorglq( n - 1, n - 1, n - 1, &A(1,1), lda, tau, work, 
                        lwork, iinfo );
            }
        }
    }
    return;
  
    //     End of DORGBR
  
#undef  A
} // end of function 

