/*
 * C++ implementation of lapack routine dorm2r
 *
 * $Id: dorm2r.cpp,v 1.5 1993/04/06 20:41:44 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:53
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dorm2r.cpp,v $
 * Revision 1.5  1993/04/06 20:41:44  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:16:19  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:08:06  alv
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

RWLAPKDECL void /*FUNCTION*/ dorm2r(const char &side, const char &trans, const long &m, const long &n, 
                                    const long &k, double *a, const long &lda, double tau[], double *c, const long &ldc, 
                                    double work[], long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
#define C(I_,J_)  (*(c+(I_)*(ldc)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
// end of PARAMETER translations

    int left, notran;
    long _do0, _do1, i1, i2, i3, i_, ic, jc, mi, ni, nq;
    long i=0;
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
  
    //  DORM2R overwrites the general real m by n matrix C with
  
    //        Q * C  if SIDE = 'L' and TRANS = 'N', or
  
    //        Q'* C  if SIDE = 'L' and TRANS = 'T', or
  
    //        C * Q  if SIDE = 'R' and TRANS = 'N', or
  
    //        C * Q' if SIDE = 'R' and TRANS = 'T',
  
    //  where Q is a real orthogonal matrix defined as the product of k
    //  elementary reflectors
  
    //        Q = H(1) H(2) . . . H(k)
  
    //  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
    //  if SIDE = 'R'.
  
    //  Arguments
    //  =========
  
    //  SIDE    (input) CHARACTER*1
    //          = 'L': apply Q or Q' from the Left
    //          = 'R': apply Q or Q' from the Right
  
    //  TRANS   (input) CHARACTER*1
    //          = 'N': apply Q  (No transpose)
    //          = 'T': apply Q' (Transpose)
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix C. M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix C. N >= 0.
  
    //  K       (input) INTEGER
    //          The number of elementary reflectors whose product defines
    //          the matrix Q.
    //          If SIDE = 'L', M >= K >= 0;
    //          if SIDE = 'R', N >= K >= 0.
  
    //  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
    //          The i-th column must contain the vector which defines the
    //          elementary reflector H(i), for i = 1,2,...,k, as returned by
    //          DGEQRF in the first k columns of its array argument A.
    //          A is modified by the routine but restored on exit.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.
    //          If SIDE = 'L', LDA >= max(1,M);
    //          if SIDE = 'R', LDA >= max(1,N).
  
    //  TAU     (input) DOUBLE PRECISION array, dimension (K)
    //          TAU(i) must contain the scalar factor of the elementary
    //          reflector H(i), as returned by DGEQRF.
  
    //  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
    //          On entry, the m by n matrix C.
    //          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
  
    //  LDC     (input) INTEGER
    //          The leading dimension of the array C. LDC >= max(1,M).
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension
    //                                   (N) if SIDE = 'L',
    //                                   (M) if SIDE = 'R'
  
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
    left = lsame( side, 'L' );
    notran = lsame( trans, 'N' );
  
    //     NQ is the order of Q
  
    if( left ) { 
        nq = m;
    }
    else { 
        nq = n;
    }
    if( !left && !lsame( side, 'R' ) ) { 
        info = -1;
    }
    else if( !notran && !lsame( trans, 'T' ) ) { 
        info = -2;
    }
    else if( m < 0 ) { 
        info = -3;
    }
    else if( n < 0 ) { 
        info = -4;
    }
    else if( k < 0 || k > nq ) { 
        info = -5;
    }
    else if( lda < max( 1, nq ) ) { 
        info = -7;
    }
    else if( ldc < max( 1, m ) ) { 
        info = -10;
    }
    if( info != 0 ) { 
        xerbla( "DORM2R", -info );
        return;
    }
  
    //     Quick return if possible
  
    if( (m == 0 || n == 0) || k == 0 ) 
        return;
  
    if( (left && !notran) || (!left && notran) ) { 
        i1 = 1;
        i2 = k;
        i3 = 1;
    }
    else { 
        i1 = k;
        i2 = 1;
        i3 = -1;
    }
  
    if( left ) { 
        ni = n;
        jc = 1;
    }
    else { 
        mi = m;
        ic = 1;
    }
  
    for( i = i1, i_ = i - 1, _do0=docnt(i,i2,_do1 = i3); _do0 > 0; i += _do1, i_ += _do1, _do0-- ) { 
        if( left ) { 
      
            //           H(i) is applied to C(i:m,1:n)
      
            mi = m - i + 1;
            ic = i;
        }
        else { 
      
            //           H(i) is applied to C(1:m,i:n)
      
            ni = n - i + 1;
            jc = i;
        }
    
        //        Apply H(i)
    
        aii = A(i_,i_);
        A(i_,i_) = ONE;
        dlarf( side, mi, ni, &A(i_,i_), 1, tau[i_], &C(jc - 1,ic - 1), 
               ldc, work );
        A(i_,i_) = aii;
    }
    return;
  
    //     End of DORM2R
  
#undef  C
#undef  A
} // end of function 

