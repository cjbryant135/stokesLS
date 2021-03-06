/*
 * C++ implementation of lapack routine dlarfb
 *
 * $Id: dlarfb.cpp,v 1.5 1993/04/06 20:41:15 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:35:53
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlarfb.cpp,v $
 * Revision 1.5  1993/04/06 20:41:15  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  18:41:23  alv
 * now passes chars explicitly, rather than indirection of a string, to shut up SUN warnings
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:42  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:35  alv
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

RWLAPKDECL void /*FUNCTION*/ dlarfb(const char &side, const char &trans, const char &direct, const char &storev, 
                                    const long &m, const long &n, const long &k, double *v, const long &ldv, 
                                    double *t, const long &ldt, double *c, const long &ldc, double *work, 
                                    const long &ldwork)
{
#define V(I_,J_)  (*(v+(I_)*(ldv)+(J_)))
#define T(I_,J_)  (*(t+(I_)*(ldt)+(J_)))
#define C(I_,J_)  (*(c+(I_)*(ldc)+(J_)))
#define WORK(I_,J_) (*(work+(I_)*(ldwork)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
// end of PARAMETER translations

    char transt;
    long _do0, _do1, _do10, _do11, _do12, _do13, _do14, _do15, 
        _do16, _do17, _do18, _do19, _do2, _do20, _do21, _do22, _do23, 
        _do3, _do4, _do5, _do6, _do7, _do8, _do9, i, i_, j, j_;

  
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
  
    //  DLARFB applies a real block reflector H or its transpose H' to a
    //  real m by n matrix C, from either the left or the right.
  
    //  Arguments
    //  =========
  
    //  SIDE    (input) CHARACTER*1
    //          = 'L': apply H or H' from the Left
    //          = 'R': apply H or H' from the Right
  
    //  TRANS   (input) CHARACTER*1
    //          = 'N': apply H (No transpose)
    //          = 'T': apply H' (Transpose)
  
    //  DIRECT  (input) CHARACTER*1
    //          Indicates how H is formed from a product of elementary
    //          reflectors
    //          = 'F': H = H(1) H(2) . . . H(k) (Forward)
    //          = 'B': H = H(k) . . . H(2) H(1) (Backward)
  
    //  STOREV  (input) CHARACTER*1
    //          Indicates how the vectors which define the elementary
    //          reflectors are stored:
    //          = 'C': Columnwise
    //          = 'R': Rowwise
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix C.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix C.
  
    //  K       (input) INTEGER
    //          The order of the matrix T (= the number of elementary
    //          reflectors whose product defines the block reflector).
  
    //  V       (input) DOUBLE PRECISION array, dimension
    //                                (LDV,K) if STOREV = 'C'
    //                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
    //                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
    //          The matrix V. See further details.
  
    //  LDV     (input) INTEGER
    //          The leading dimension of the array V.
    //          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
    //          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
    //          if STOREV = 'R', LDV >= K.
  
    //  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
    //          The triangular k by k matrix T in the representation of the
    //          block reflector.
  
    //  LDT     (input) INTEGER
    //          The leading dimension of the array T. LDT >= K.
  
    //  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
    //          On entry, the m by n matrix C.
    //          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
  
    //  LDC     (input) INTEGER
    //          The leading dimension of the array C. LDA >= max(1,M).
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
  
    //  LDWORK  (input) INTEGER
    //          The leading dimension of the array WORK.
    //          If SIDE = 'L', LDWORK >= max(1,N);
    //          if SIDE = 'R', LDWORK >= max(1,M).
  
    //  =====================================================================
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
  
    //     Quick return if possible
  
    if( m <= 0 || n <= 0 ) 
        return;
  
    if( lsame( trans, 'N' ) ) { 
        transt = 'T';
    }
    else { 
        transt = 'N';
    }
  
    if( lsame( storev, 'C' ) ) { 
    
        if( lsame( direct, 'F' ) ) { 
      
            //           Let  V =  ( V1 )    (first K rows)
            //                     ( V2 )
            //           where  V1  is unit lower triangular.
      
            if( lsame( side, 'L' ) ) { 
        
                //              Form  H * C  or  H' * C  where  C = ( C1 )
                //                                                  ( C2 )
        
                //              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
        
                //              W := C1'
        
                for( j = 1, j_ = j - 1, _do0 = k; j <= _do0; j++, j_++ ) { 
                    dcopy( n, &C(0,j_), ldc, &WORK(j_,0), 1 );
                }
        
                //              W := W * V1
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'N'/* No transpose */
                       , 'U'/* Unit */, n, k, ONE, v, ldv, work, ldwork );
                if( m > k ) { 
          
                    //                 W := W + C2'*V2
          
                    dgemm( 'T'/* Transpose */, 'N'/* No transpose */
                           , n, k, m - k, ONE, &C(0,k), ldc, &V(0,k), ldv, 
                           ONE, work, ldwork );
                }
        
                //              W := W * T'  or  W * T
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, transt, 'N'/* Non-unit */
                       , n, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - V * W'
        
                if( m > k ) { 
          
                    //                 C2 := C2 - V2 * W'
          
                    dgemm( 'N'/* No transpose */, 'T'/* Transpose */
                           , m - k, n, k, -ONE, &V(0,k), ldv, work, ldwork, 
                           ONE, &C(0,k), ldc );
                }
        
                //              W := W * V1'
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'T'/* Transpose */
                       , 'U'/* Unit */, n, k, ONE, v, ldv, work, ldwork );
        
                //              C1 := C1 - W'
        
                for( j = 1, j_ = j - 1, _do1 = k; j <= _do1; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do2 = n; i <= _do2; i++, i_++ ) { 
                        C(i_,j_) = C(i_,j_) - WORK(j_,i_);
                    }
                }
        
            }
            else if( lsame( side, 'R' ) ) { 
        
                //              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        
                //              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
        
                //              W := C1
        
                for( j = 1, j_ = j - 1, _do3 = k; j <= _do3; j++, j_++ ) { 
                    dcopy( m, &C(j_,0), 1, &WORK(j_,0), 1 );
                }
        
                //              W := W * V1
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'N'/* No transpose */
                       , 'U'/* Unit */, m, k, ONE, v, ldv, work, ldwork );
                if( n > k ) { 
          
                    //                 W := W + C2 * V2
          
                    dgemm( 'N'/* No transpose */, 'N'/* No transpose */
                           , m, k, n - k, ONE, &C(k,0), ldc, &V(0,k), ldv, 
                           ONE, work, ldwork );
                }
        
                //              W := W * T  or  W * T'
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, trans, 'N'/* Non-unit */
                       , m, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - W * V'
        
                if( n > k ) { 
          
                    //                 C2 := C2 - W * V2'
          
                    dgemm( 'N'/* No transpose */, 'T'/* Transpose */
                           , m, n - k, k, -ONE, work, ldwork, &V(0,k), ldv, 
                           ONE, &C(k,0), ldc );
                }
        
                //              W := W * V1'
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'T'/* Transpose */
                       , 'U'/* Unit */, m, k, ONE, v, ldv, work, ldwork );
        
                //              C1 := C1 - W
        
                for( j = 1, j_ = j - 1, _do4 = k; j <= _do4; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do5 = m; i <= _do5; i++, i_++ ) { 
                        C(j_,i_) = C(j_,i_) - WORK(j_,i_);
                    }
                }
            }
      
        }
        else { 
      
            //           Let  V =  ( V1 )
            //                     ( V2 )    (last K rows)
            //           where  V2  is unit upper triangular.
      
            if( lsame( side, 'L' ) ) { 
        
                //              Form  H * C  or  H' * C  where  C = ( C1 )
                //                                                  ( C2 )
        
                //              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
        
                //              W := C2'
        
                for( j = 1, j_ = j - 1, _do6 = k; j <= _do6; j++, j_++ ) { 
                    dcopy( n, &C(0,m - k + j_), ldc, &WORK(j_,0), 
                           1 );
                }
        
                //              W := W * V2
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'N'/* No transpose */
                       , 'U'/* Unit */, n, k, ONE, &V(0,m - k), ldv, work, 
                       ldwork );
                if( m > k ) { 
          
                    //                 W := W + C1'*V1
          
                    dgemm( 'T'/* Transpose */, 'N'/* No transpose */
                           , n, k, m - k, ONE, c, ldc, v, ldv, ONE, work, 
                           ldwork );
                }
        
                //              W := W * T'  or  W * T
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, transt, 'N'/* Non-unit */
                       , n, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - V * W'
        
                if( m > k ) { 
          
                    //                 C1 := C1 - V1 * W'
          
                    dgemm( 'N'/* No transpose */, 'T'/* Transpose */
                           , m - k, n, k, -ONE, v, ldv, work, ldwork, ONE, 
                           c, ldc );
                }
        
                //              W := W * V2'
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'T'/* Transpose */
                       , 'U'/* Unit */, n, k, ONE, &V(0,m - k), ldv, work, 
                       ldwork );
        
                //              C2 := C2 - W'
        
                for( j = 1, j_ = j - 1, _do7 = k; j <= _do7; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do8 = n; i <= _do8; i++, i_++ ) { 
                        C(i_,m - k + j_) = C(i_,m - k + j_) - WORK(j_,i_);
                    }
                }
        
            }
            else if( lsame( side, 'R' ) ) { 
        
                //              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        
                //              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
        
                //              W := C2
        
                for( j = 1, j_ = j - 1, _do9 = k; j <= _do9; j++, j_++ ) { 
                    dcopy( m, &C(n - k + j_,0), 1, &WORK(j_,0), 1 );
                }
        
                //              W := W * V2
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'N'/* No transpose */
                       , 'U'/* Unit */, m, k, ONE, &V(0,n - k), ldv, work, 
                       ldwork );
                if( n > k ) { 
          
                    //                 W := W + C1 * V1
          
                    dgemm( 'N'/* No transpose */, 'N'/* No transpose */
                           , m, k, n - k, ONE, c, ldc, v, ldv, ONE, work, 
                           ldwork );
                }
        
                //              W := W * T  or  W * T'
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, trans, 'N'/* Non-unit */
                       , m, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - W * V'
        
                if( n > k ) { 
          
                    //                 C1 := C1 - W * V1'
          
                    dgemm( 'N'/* No transpose */, 'T'/* Transpose */
                           , m, n - k, k, -ONE, work, ldwork, v, ldv, ONE, 
                           c, ldc );
                }
        
                //              W := W * V2'
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'T'/* Transpose */
                       , 'U'/* Unit */, m, k, ONE, &V(0,n - k), ldv, work, 
                       ldwork );
        
                //              C2 := C2 - W
        
                for( j = 1, j_ = j - 1, _do10 = k; j <= _do10; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do11 = m; i <= _do11; i++, i_++ ) { 
                        C(n - k + j_,i_) = C(n - k + j_,i_) - WORK(j_,i_);
                    }
                }
            }
        }
    
    }
    else if( lsame( storev, 'R' ) ) { 
    
        if( lsame( direct, 'F' ) ) { 
      
            //           Let  V =  ( V1  V2 )    (V1: first K columns)
            //           where  V1  is unit upper triangular.
      
            if( lsame( side, 'L' ) ) { 
        
                //              Form  H * C  or  H' * C  where  C = ( C1 )
                //                                                  ( C2 )
        
                //              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
        
                //              W := C1'
        
                for( j = 1, j_ = j - 1, _do12 = k; j <= _do12; j++, j_++ ) { 
                    dcopy( n, &C(0,j_), ldc, &WORK(j_,0), 1 );
                }
        
                //              W := W * V1'
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'T'/* Transpose */
                       , 'U'/* Unit */, n, k, ONE, v, ldv, work, ldwork );
                if( m > k ) { 
          
                    //                 W := W + C2'*V2'
          
                    dgemm( 'T'/* Transpose */, 'T'/* Transpose */, n, 
                           k, m - k, ONE, &C(0,k), ldc, &V(k,0), ldv, ONE, 
                           work, ldwork );
                }
        
                //              W := W * T'  or  W * T
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, transt, 'N'/* Non-unit */
                       , n, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - V' * W'
        
                if( m > k ) { 
          
                    //                 C2 := C2 - V2' * W'
          
                    dgemm( 'T'/* Transpose */, 'T'/* Transpose */, m - 
                           k, n, k, -ONE, &V(k,0), ldv, work, ldwork, ONE, 
                           &C(0,k), ldc );
                }
        
                //              W := W * V1
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'N'/* No transpose */
                       , 'U'/* Unit */, n, k, ONE, v, ldv, work, ldwork );
        
                //              C1 := C1 - W'
        
                for( j = 1, j_ = j - 1, _do13 = k; j <= _do13; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do14 = n; i <= _do14; i++, i_++ ) { 
                        C(i_,j_) = C(i_,j_) - WORK(j_,i_);
                    }
                }
        
            }
            else if( lsame( side, 'R' ) ) { 
        
                //              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        
                //              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
        
                //              W := C1
        
                for( j = 1, j_ = j - 1, _do15 = k; j <= _do15; j++, j_++ ) { 
                    dcopy( m, &C(j_,0), 1, &WORK(j_,0), 1 );
                }
        
                //              W := W * V1'
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'T'/* Transpose */
                       , 'U'/* Unit */, m, k, ONE, v, ldv, work, ldwork );
                if( n > k ) { 
          
                    //                 W := W + C2 * V2'
          
                    dgemm( 'N'/* No transpose */, 'T'/* Transpose */
                           , m, k, n - k, ONE, &C(k,0), ldc, &V(k,0), ldv, 
                           ONE, work, ldwork );
                }
        
                //              W := W * T  or  W * T'
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, trans, 'N'/* Non-unit */
                       , m, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - W * V
        
                if( n > k ) { 
          
                    //                 C2 := C2 - W * V2
          
                    dgemm( 'N'/* No transpose */, 'N'/* No transpose */
                           , m, n - k, k, -ONE, work, ldwork, &V(k,0), ldv, 
                           ONE, &C(k,0), ldc );
                }
        
                //              W := W * V1
        
                dtrmm( 'R'/* Right */, 'U'/* Upper */, 'N'/* No transpose */
                       , 'U'/* Unit */, m, k, ONE, v, ldv, work, ldwork );
        
                //              C1 := C1 - W
        
                for( j = 1, j_ = j - 1, _do16 = k; j <= _do16; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do17 = m; i <= _do17; i++, i_++ ) { 
                        C(j_,i_) = C(j_,i_) - WORK(j_,i_);
                    }
                }
        
            }
      
        }
        else { 
      
            //           Let  V =  ( V1  V2 )    (V2: last K columns)
            //           where  V2  is unit lower triangular.
      
            if( lsame( side, 'L' ) ) { 
        
                //              Form  H * C  or  H' * C  where  C = ( C1 )
                //                                                  ( C2 )
        
                //              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
        
                //              W := C2'
        
                for( j = 1, j_ = j - 1, _do18 = k; j <= _do18; j++, j_++ ) { 
                    dcopy( n, &C(0,m - k + j_), ldc, &WORK(j_,0), 
                           1 );
                }
        
                //              W := W * V2'
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'T'/* Transpose */
                       , 'U'/* Unit */, n, k, ONE, &V(m - k,0), ldv, work, 
                       ldwork );
                if( m > k ) { 
          
                    //                 W := W + C1'*V1'
          
                    dgemm( 'T'/* Transpose */, 'T'/* Transpose */, n, 
                           k, m - k, ONE, c, ldc, v, ldv, ONE, work, ldwork );
                }
        
                //              W := W * T'  or  W * T
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, transt, 'N'/* Non-unit */
                       , n, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - V' * W'
        
                if( m > k ) { 
          
                    //                 C1 := C1 - V1' * W'
          
                    dgemm( 'T'/* Transpose */, 'T'/* Transpose */, m - 
                           k, n, k, -ONE, v, ldv, work, ldwork, ONE, c, 
                           ldc );
                }
        
                //              W := W * V2
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'N'/* No transpose */
                       , 'U'/* Unit */, n, k, ONE, &V(m - k,0), ldv, work, 
                       ldwork );
        
                //              C2 := C2 - W'
        
                for( j = 1, j_ = j - 1, _do19 = k; j <= _do19; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do20 = n; i <= _do20; i++, i_++ ) { 
                        C(i_,m - k + j_) = C(i_,m - k + j_) - WORK(j_,i_);
                    }
                }
        
            }
            else if( lsame( side, 'R' ) ) { 
        
                //              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        
                //              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
        
                //              W := C2
        
                for( j = 1, j_ = j - 1, _do21 = k; j <= _do21; j++, j_++ ) { 
                    dcopy( m, &C(n - k + j_,0), 1, &WORK(j_,0), 1 );
                }
        
                //              W := W * V2'
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'T'/* Transpose */
                       , 'U'/* Unit */, m, k, ONE, &V(n - k,0), ldv, work, 
                       ldwork );
                if( n > k ) { 
          
                    //                 W := W + C1 * V1'
          
                    dgemm( 'N'/* No transpose */, 'T'/* Transpose */
                           , m, k, n - k, ONE, c, ldc, v, ldv, ONE, work, 
                           ldwork );
                }
        
                //              W := W * T  or  W * T'
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, trans, 'N'/* Non-unit */
                       , m, k, ONE, t, ldt, work, ldwork );
        
                //              C := C - W * V
        
                if( n > k ) { 
          
                    //                 C1 := C1 - W * V1
          
                    dgemm( 'N'/* No transpose */, 'N'/* No transpose */
                           , m, n - k, k, -ONE, work, ldwork, v, ldv, ONE, 
                           c, ldc );
                }
        
                //              W := W * V2
        
                dtrmm( 'R'/* Right */, 'L'/* Lower */, 'N'/* No transpose */
                       , 'U'/* Unit */, m, k, ONE, &V(n - k,0), ldv, work, 
                       ldwork );
        
                //              C1 := C1 - W
        
                for( j = 1, j_ = j - 1, _do22 = k; j <= _do22; j++, j_++ ) { 
                    for( i = 1, i_ = i - 1, _do23 = m; i <= _do23; i++, i_++ ) { 
                        C(n - k + j_,i_) = C(n - k + j_,i_) - WORK(j_,i_);
                    }
                }
        
            }
      
        }
    }
  
    return;
  
    //     End of DLARFB
  
#undef  WORK
#undef  V
#undef  T
#undef  C
} // end of function 

