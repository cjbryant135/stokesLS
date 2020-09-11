/*
 * Default C++ implementation of dger
 * For optimum performance, use a machine specific bla library
 *
 * $Id: dger.cpp,v 1.4 1993/03/19 16:26:58 alv Exp $
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
 *
 * Rewritten in C++ for efficiency using unrolled loops.
 *
 * $Log: dger.cpp,v $
 * Revision 1.4  1993/03/19 16:26:58  alv
 * added RWBLADECL linkage specification
 *
 * Revision 1.3  1993/03/05  23:07:15  alv
 * changed ref parms to const ref
 *
 * Revision 1.2  1993/03/05  22:41:24  alv
 * rewritten in efficient C++ with unrolled loops
 *
 * Revision 1.1  1993/03/03  16:04:55  alv
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

RWLAPKDECL void dger(const long &m, const long &n, const double &alpha, double x[], 
                     const long &incx, double y[], const long &incy, double *a, const long &lda)
{
    long N=n;
    long info;

    //     .. Scalar Arguments ..
    //     .. Array Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  DGER   performs the rank 1 operation
  
    //     A := alpha*x*y' + A,
  
    //  where alpha is a scalar, x is an m element vector, y is an n element
    //  vector and A is an m by n matrix.
  
    //  Parameters
    //  ==========
  
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
  
    //  X      - DOUBLE PRECISION array of dimension at least
    //           ( 1 + ( m - 1 )*abs( INCX ) ).
    //           Before entry, the incremented array X must contain the m
    //           element vector x.
    //           Unchanged on exit.
  
    //  INCX   - INTEGER.
    //           On entry, INCX specifies the increment for the elements of
    //           X. INCX must not be zero.
    //           Unchanged on exit.
  
    //  Y      - DOUBLE PRECISION array of dimension at least
    //           ( 1 + ( n - 1 )*abs( INCY ) ).
    //           Before entry, the incremented array Y must contain the n
    //           element vector y.
    //           Unchanged on exit.
  
    //  INCY   - INTEGER.
    //           On entry, INCY specifies the increment for the elements of
    //           Y. INCY must not be zero.
    //           Unchanged on exit.
  
    //  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
    //           Before entry, the leading m by n part of the array A must
    //           contain the matrix of coefficients. On exit, A is
    //           overwritten by the updated matrix.
  
    //  LDA    - INTEGER.
    //           On entry, LDA specifies the first dimension of A as declared
    //           in the calling (sub) program. LDA must be at least
    //           max( 1, m ).
    //           Unchanged on exit.
  
  
    //  Level 2 Blas routine.
  
    //  -- Written on 22-October-1986.
    //     Jack Dongarra, Argonne National Lab.
    //     Jeremy Du Croz, Nag Central Office.
    //     Sven Hammarling, Nag Central Office.
    //     Richard Hanson, Sandia National Labs.
  
  
    //     .. Parameters ..
    //     .. Local Scalars ..
    //     .. External Subroutines ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test the input parameters.
  
    info = 0;
    if( m < 0 ) { 
        info = 1;
    }
    else if( n < 0 ) { 
        info = 2;
    }
    else if( incx == 0 ) { 
        info = 5;
    }
    else if( incy == 0 ) { 
        info = 7;
    }
    else if( lda < max( 1, m ) ) { 
        info = 9;
    }
    if( info != 0 ) { 
        xerbla( "DGER  ", info );
        return;
    }
  
    //     Quick return if possible.
  
    if( ((m == 0) || (n == 0)) || (alpha == 0) ) 
        return;
  
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through A.
  
    if( incy < 0 ) { 
        y += (n-1)*incy;
    }
    if ( incx < 0 ) {
        x += (m-1)*incx;
    }
    if (incx==1) {
        while(N--) {
            if (*y != 0) {
                double temp = alpha*(*y);
                const double *xptr = x;
                double *aptr = a;
                int i;
                for(i=m; i>=5; i-=5) {
                    aptr[0] += temp*xptr[0];
                    aptr[1] += temp*xptr[1];
                    aptr[2] += temp*xptr[2];
                    aptr[3] += temp*xptr[3];
                    aptr[4] += temp*xptr[4];
                    aptr += 5;
                    xptr += 5;
                }
                while (i--) {
                    *aptr++ += temp*(*xptr++);
                }
            }
            a += lda;
            y += incy;
        }
    } else {
        while(N--) {
            if (*y != 0) {
                double temp = alpha*(*y);
                const double *xptr = x;
                double *aptr = a;
                for(register int i=m; i--;) {
                    *aptr++ += temp*(*xptr);
                    xptr += incx;
                }
            }
            a += lda;
            y += incy;
        }
    }
}

