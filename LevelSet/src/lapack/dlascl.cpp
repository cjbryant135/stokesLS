/*
 * C++ implementation of lapack routine dlascl
 *
 * $Id: dlascl.cpp,v 1.5 1993/04/06 20:41:26 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:09
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlascl.cpp,v $
 * Revision 1.5  1993/04/06 20:41:26  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:15:52  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:42  alv
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

RWLAPKDECL void /*FUNCTION*/ dlascl(const char &type, const long &kl, const long &ku, const double &cfrom, 
                                    const double &cto, const long &m, const long &n, double *a, const long &lda, 
                                    long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ZERO = 0.0e0;
    const double ONE = 1.0e0;
// end of PARAMETER translations

    int done;
    long _do0, _do1, _do10, _do11, _do12, _do13, _do2, _do3, _do4, 
        _do5, _do6, _do7, _do8, _do9, i, i_, itype, j, j_, k1, k2, k3, 
        k4;
    double bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;

  
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
  
    //  DLASCL multiplies the M by N real matrix A by the real scalar
    //  CTO/CFROM.  This is done without over/underflow as long as the final
    //  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
    //  A may be full, upper triangular, lower triangular, upper Hessenberg,
    //  or banded.
  
    //  Arguments
    //  =========
  
    //  TYPE    (input) CHARACTER*1
    //          TYPE indices the storage type of the input matrix.
    //          = 'G':  A is a full matrix.
    //          = 'L':  A is a lower triangular matrix.
    //          = 'U':  A is an upper triangular matrix.
    //          = 'H':  A is an upper Hessenberg matrix.
    //          = 'B':  A is a symmetric band matrix with lower bandwidth KL
    //                  and upper bandwidth KU and with the only the lower
    //                  half stored.
    //          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
    //                  and upper bandwidth KU and with the only the upper
    //                  half stored.
    //          = 'Z':  A is a band matrix with lower bandwidth KL and upper
    //                  bandwidth KU.
  
    //  KL      (input) INTEGER
    //          The lower bandwidth of A.  Referenced only if TYPE = 'B',
    //          'Q' or 'Z'.
  
    //  KU      (input) INTEGER
    //          The upper bandwidth of A.  Referenced only if TYPE = 'B',
    //          'Q' or 'Z'.
  
    //  CFROM   (input) DOUBLE PRECISION
    //  CTO     (input) DOUBLE PRECISION
    //          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
    //          without over/underflow if the final result CTO*A(I,J)/CFROM
    //          can be represented without over/underflow.  CFROM must be
    //          nonzero.
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix A.  M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix A.  N >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
    //          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
    //          storage type.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
    //  INFO    (output) INTEGER
    //          0  - successful exit
    //          <0 - if INFO = -i, the i-th argument had an illegal value.
  
    //  =====================================================================
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test the input arguments
  
    info = 0;
  
    if( lsame( type, 'G' ) ) { 
        itype = 0;
    }
    else if( lsame( type, 'L' ) ) { 
        itype = 1;
    }
    else if( lsame( type, 'U' ) ) { 
        itype = 2;
    }
    else if( lsame( type, 'H' ) ) { 
        itype = 3;
    }
    else if( lsame( type, 'B' ) ) { 
        itype = 4;
    }
    else if( lsame( type, 'Q' ) ) { 
        itype = 5;
    }
    else if( lsame( type, 'Z' ) ) { 
        itype = 6;
    }
    else { 
        itype = -1;
    }
  
    if( itype == -1 ) { 
        info = -1;
    }
    else if( cfrom == ZERO ) { 
        info = -4;
    }
    else if( m < 0 ) { 
        info = -6;
    }
    else if( (n < 0 || (itype == 4 && n != m)) || (itype == 5 && n != 
                                                   m) ) { 
        info = -7;
    }
    else if( itype <= 3 && lda < max( 1, m ) ) { 
        info = -9;
    }
    else if( itype >= 4 ) { 
        if( kl < 0 || kl > max( m - 1, 0 ) ) { 
            info = -2;
        }
        else if( (ku < 0 || ku > max( n - 1, 0 )) || ((itype == 4 || 
                                                       itype == 5) && kl != ku) ) { 
            info = -3;
        }
        else if( ((itype == 4 && lda < kl + 1) || (itype == 5 && lda < 
                                                   ku + 1)) || (itype == 6 && lda < 2*kl + ku + 1) ) { 
            info = -9;
        }
    }
  
    if( info != 0 ) { 
        xerbla( "DLASCL", -info );
        return;
    }
  
    //     Quick return if possible
  
    if( n == 0 || m == 0 ) 
        return;
  
    //     Get machine parameters
  
    smlnum = dlamch( 'S' );
    bignum = ONE/smlnum;
  
    cfromc = cfrom;
    ctoc = cto;
  
L_10:
    ;
    cfrom1 = cfromc*smlnum;
    cto1 = ctoc/bignum;
    if( fabs( cfrom1 ) > fabs( ctoc ) && ctoc != ZERO ) { 
        mul = smlnum;
        done = FALSE;
        cfromc = cfrom1;
    }
    else if( fabs( cto1 ) > fabs( cfromc ) ) { 
        mul = bignum;
        done = FALSE;
        ctoc = cto1;
    }
    else { 
        mul = ctoc/cfromc;
        done = TRUE;
    }
  
    if( itype == 0 ) { 
    
        //        Full matrix
    
        for( j = 1, j_ = j - 1, _do0 = n; j <= _do0; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do1 = m; i <= _do1; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
    else if( itype == 1 ) { 
    
        //        Lower triangular matrix
    
        for( j = 1, j_ = j - 1, _do2 = n; j <= _do2; j++, j_++ ) { 
            for( i = j, i_ = i - 1, _do3 = m; i <= _do3; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
    else if( itype == 2 ) { 
    
        //        Upper triangular matrix
    
        for( j = 1, j_ = j - 1, _do4 = n; j <= _do4; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do5 = min( j, m ); i <= _do5; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
    else if( itype == 3 ) { 
    
        //        Upper Hessenberg matrix
    
        for( j = 1, j_ = j - 1, _do6 = n; j <= _do6; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do7 = min( j + 1, m ); i <= _do7; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
    else if( itype == 4 ) { 
    
        //        Lower half of a symmetric band matrix
    
        k3 = kl + 1;
        k4 = n + 1;
        for( j = 1, j_ = j - 1, _do8 = n; j <= _do8; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do9 = min( k3, k4 - j ); i <= _do9; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
    else if( itype == 5 ) { 
    
        //        Upper half of a symmetric band matrix
    
        k1 = ku + 2;
        k3 = ku + 1;
        for( j = 1, j_ = j - 1, _do10 = n; j <= _do10; j++, j_++ ) { 
            for( i = max( k1 - j, 1 ), i_ = i - 1, _do11 = k3; i <= _do11; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
    else if( itype == 6 ) { 
    
        //        Band matrix
    
        k1 = kl + ku + 2;
        k2 = kl + 1;
        k3 = 2*kl + ku + 1;
        k4 = kl + ku + 1 + m;
        for( j = 1, j_ = j - 1, _do12 = n; j <= _do12; j++, j_++ ) { 
            for( i = max( k1 - j, k2 ), i_ = i - 1, _do13 = min( k3, 
                                                                 k4 - j ); i <= _do13; i++, i_++ ) { 
                A(j_,i_) = A(j_,i_)*mul;
            }
        }
    
    }
  
    if( !done ) 
        goto L_10;
  
    return;
  
    //     End of DLASCL
  
#undef  A
} // end of function 

