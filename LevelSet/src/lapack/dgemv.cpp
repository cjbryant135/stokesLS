/*
 * Default C++ implementation of dgemv
 * For optimum performance, use a machine specific bla library
 *
 * $Id: dgemv.cpp,v 1.3 1993/03/19 16:26:58 alv Exp $
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
 * Translated from the default Fortran bla using Cobalt
 * Blue's FOR_C++, and then massaged slightly to Rogue
 * Wave format.
 *
 * Translated by FOR_C++, v1.1 (P), on 02/17/93 at 14:40:06
 * FOR_C++ Options SET: alloc do=rt no=p pf=dbla,xbla s=dv str=l - prototypes
 *
 * $Log: dgemv.cpp,v $
 * Revision 1.3  1993/03/19 16:26:58  alv
 * added RWBLADECL linkage specification
 *
 * Revision 1.2  1993/03/05  23:07:14  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:04:36  alv
 * Initial revision
 *
 */

#define RW_CPPBLAS 1
#if 0
#include "rw/bla.h"
#include "rw/bla.h"
#include "rw/fortran.h" /* Fortran run time library */
#else
#include "level/lapack.h"
#endif

RWLAPKDECL void dgemv(const char &trans, const long &m, const long &n, const double &alpha, 
                      double *a, const long &lda, double x[], const long &incx, const double &beta, 
                      double y[], const long &incy)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, _do10, _do11, _do2, _do3, _do4, _do5, _do6, 
        _do7, _do8, _do9, i, i_, info, ix, iy, j, j_, jx, jy, kx, ky, 
        lenx, leny;
    double temp;

    //     .. Scalar Arguments ..
    //     .. Array Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  DGEMV  performs one of the matrix-vector operations
  
    //     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  
    //  where alpha and beta are scalars, x and y are vectors and A is an
    //  m by n matrix.
  
    //  Parameters
    //  ==========
  
    //  TRANS  - CHARACTER*1.
    //           On entry, TRANS specifies the operation to be performed as
    //           follows:
  
    //              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
  
    //              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
  
    //              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
  
    //           Unchanged on exit.
  
    //  M      - INTEGER.
    //           On entry, M specifies the number of rows of the matrix A.
    //           M must be at least zero.
    //           Unchanged on exit.
  
    //  N      - INTEGER.
    //           On entry, N specifies the number of columns of the matrix A.
    //           N must be at least zero.
    //           Unchanged on exit.
  
    //  ALPHA  - DOUBLE PRECISION.
    //           On entry, ALPHA specifies the scalar alpha.
    //           Unchanged on exit.
  
    //  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
    //           Before entry, the leading m by n part of the array A must
    //           contain the matrix of coefficients.
    //           Unchanged on exit.
  
    //  LDA    - INTEGER.
    //           On entry, LDA specifies the first dimension of A as declared
    //           in the calling (sub) program. LDA must be at least
    //           max( 1, m ).
    //           Unchanged on exit.
  
    //  X      - DOUBLE PRECISION array of DIMENSION at least
    //           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
    //           and at least
    //           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
    //           Before entry, the incremented array X must contain the
    //           vector x.
    //           Unchanged on exit.
  
    //  INCX   - INTEGER.
    //           On entry, INCX specifies the increment for the elements of
    //           X. INCX must not be zero.
    //           Unchanged on exit.
  
    //  BETA   - DOUBLE PRECISION.
    //           On entry, BETA specifies the scalar beta. When BETA is
    //           supplied as zero then Y need not be set on input.
    //           Unchanged on exit.
  
    //  Y      - DOUBLE PRECISION array of DIMENSION at least
    //           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
    //           and at least
    //           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
    //           Before entry with BETA non-zero, the incremented array Y
    //           must contain the vector y. On exit, Y is overwritten by the
    //           updated vector y.
  
    //  INCY   - INTEGER.
    //           On entry, INCY specifies the increment for the elements of
    //           Y. INCY must not be zero.
    //           Unchanged on exit.
  
  
    //  Level 2 Blas routine.
  
    //  -- Written on 22-October-1986.
    //     Jack Dongarra, Argonne National Lab.
    //     Jeremy Du Croz, Nag Central Office.
    //     Sven Hammarling, Nag Central Office.
    //     Richard Hanson, Sandia National Labs.
  
  
    //     .. Parameters ..
    //     .. Local Scalars ..
    //     .. External Functions ..
    //     .. External Subroutines ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test the input parameters.
  
    info = 0;
    if( (!lsame( trans, 'N' ) && !lsame( trans, 'T' )) && !lsame( trans, 
                                                                  'C' ) ) { 
        info = 1;
    }
    else if( m < 0 ) { 
        info = 2;
    }
    else if( n < 0 ) { 
        info = 3;
    }
    else if( lda < max( 1, m ) ) { 
        info = 6;
    }
    else if( incx == 0 ) { 
        info = 8;
    }
    else if( incy == 0 ) { 
        info = 11;
    }
    if( info != 0 ) { 
        xerbla( "DGEMV ", info );
        return;
    }
  
    //     Quick return if possible.
  
    if( ((m == 0) || (n == 0)) || ((alpha == ZERO) && (beta == ONE)) ) 
        return;
  
    //     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    //     up the start points in  X  and  Y.
  
    if( lsame( trans, 'N' ) ) { 
        lenx = n;
        leny = m;
    }
    else { 
        lenx = m;
        leny = n;
    }
    if( incx > 0 ) { 
        kx = 1;
    }
    else { 
        kx = 1 - (lenx - 1)*incx;
    }
    if( incy > 0 ) { 
        ky = 1;
    }
    else { 
        ky = 1 - (leny - 1)*incy;
    }
  
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through A.
  
    //     First form  y := beta*y.
  
    if( beta != ONE ) { 
        if( incy == 1 ) { 
            if( beta == ZERO ) { 
                for( i = 1, i_ = i - 1, _do0 = leny; i <= _do0; i++, i_++ ) { 
                    y[i_] = ZERO;
                }
            }
            else { 
                for( i = 1, i_ = i - 1, _do1 = leny; i <= _do1; i++, i_++ ) { 
                    y[i_] = beta*y[i_];
                }
            }
        }
        else { 
            iy = ky;
            if( beta == ZERO ) { 
                for( i = 1, i_ = i - 1, _do2 = leny; i <= _do2; i++, i_++ ) { 
                    y[iy - 1] = ZERO;
                    iy = iy + incy;
                }
            }
            else { 
                for( i = 1, i_ = i - 1, _do3 = leny; i <= _do3; i++, i_++ ) { 
                    y[iy - 1] = beta*y[iy - 1];
                    iy = iy + incy;
                }
            }
        }
    }
    if( alpha == ZERO ) 
        return;
    if( lsame( trans, 'N' ) ) { 
    
        //        Form  y := alpha*A*x + y.
    
        jx = kx;
        if( incy == 1 ) { 
            for( j = 1, j_ = j - 1, _do4 = n; j <= _do4; j++, j_++ ) { 
                if( x[jx - 1] != ZERO ) { 
                    temp = alpha*x[jx - 1];
                    for( i = 1, i_ = i - 1, _do5 = m; i <= _do5; i++, i_++ ) { 
                        y[i_] = y[i_] + temp*A(j_,i_);
                    }
                }
                jx = jx + incx;
            }
        }
        else { 
            for( j = 1, j_ = j - 1, _do6 = n; j <= _do6; j++, j_++ ) { 
                if( x[jx - 1] != ZERO ) { 
                    temp = alpha*x[jx - 1];
                    iy = ky;
                    for( i = 1, i_ = i - 1, _do7 = m; i <= _do7; i++, i_++ ) { 
                        y[iy - 1] = y[iy - 1] + temp*A(j_,i_);
                        iy = iy + incy;
                    }
                }
                jx = jx + incx;
            }
        }
    }
    else { 
    
        //        Form  y := alpha*A'*x + y.
    
        jy = ky;
        if( incx == 1 ) { 
            for( j = 1, j_ = j - 1, _do8 = n; j <= _do8; j++, j_++ ) { 
                temp = ZERO;
                for( i = 1, i_ = i - 1, _do9 = m; i <= _do9; i++, i_++ ) { 
                    temp = temp + A(j_,i_)*x[i_];
                }
                y[jy - 1] = y[jy - 1] + alpha*temp;
                jy = jy + incy;
            }
        }
        else { 
            for( j = 1, j_ = j - 1, _do10 = n; j <= _do10; j++, j_++ ) { 
                temp = ZERO;
                ix = kx;
                for( i = 1, i_ = i - 1, _do11 = m; i <= _do11; i++, i_++ ) { 
                    temp = temp + A(j_,i_)*x[ix - 1];
                    ix = ix + incx;
                }
                y[jy - 1] = y[jy - 1] + alpha*temp;
                jy = jy + incy;
            }
        }
    }
  
    return;
  
    //     End of DGEMV .
  
#undef  A
} // end of function 

