/*
 * C++ implementation of lapack routine dlabrd
 *
 * $Id: dlabrd.cpp,v 1.5 1993/04/06 20:40:52 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:34:56
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlabrd.cpp,v $
 * Revision 1.5  1993/04/06 20:40:52  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  18:41:23  alv
 * now passes chars explicitly, rather than indirection of a string, to shut up SUN warnings
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:06  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:03  alv
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

RWLAPKDECL void /*FUNCTION*/ dlabrd(const long &m, const long &n, const long &nb, double *a, 
                                    const long &lda, double d[], double e[], double tauq[], double taup[], 
                                    double *x, const long &ldx, double *y, const long &ldy)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
#define X(I_,J_)  (*(x+(I_)*(ldx)+(J_)))
#define Y(I_,J_)  (*(y+(I_)*(ldy)+(J_)))
// PARAMETER translations
    const double ZERO = 0.0e0;
    const double ONE = 1.0e0;
// end of PARAMETER translations

    long _do0, _do1, i, i_;

  
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
  
    //  DLABRD reduces the first NB rows and columns of a real general
    //  m by n matrix A to upper or lower bidiagonal form by an orthogonal
    //  transformation Q' * A * P, and returns the matrices X and Y which
    //  are needed to apply the transformation to the unreduced part of A.
  
    //  If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
    //  bidiagonal form.
  
    //  This is an auxiliary routine called by DGEBRD
  
    //  Arguments
    //  =========
  
    //  M       (input) INTEGER
    //          The number of rows in the matrix A.
  
    //  N       (input) INTEGER
    //          The number of columns in the matrix A.
  
    //  NB      (input) INTEGER
    //          The number of leading rows and columns of A to be reduced.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the m by n general matrix to be reduced.
    //          On exit, the first NB rows and columns of the matrix are
    //          overwritten; the rest of the array is unchanged.
    //          If m >= n, elements on and below the diagonal in the first NB
    //            columns, with the array TAUQ, represent the orthogonal
    //            matrix Q as a product of elementary reflectors; and
    //            elements above the diagonal in the first NB rows, with the
    //            array TAUP, represent the orthogonal matrix P as a product
    //            of elementary reflectors.
    //          If m < n, elements below the diagonal in the first NB
    //            columns, with the array TAUQ, represent the orthogonal
    //            matrix Q as a product of elementary reflectors, and
    //            elements on and above the diagonal in the first NB rows,
    //            with the array TAUP, represent the orthogonal matrix P as
    //            a product of elementary reflectors.
    //          See Further Details.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
    //  D       (output) DOUBLE PRECISION array, dimension (NB)
    //          The diagonal elements of the first NB rows and columns of
    //          the reduced matrix.  D(i) = A(i,i).
  
    //  E       (output) DOUBLE PRECISION array, dimension (NB)
    //          The off-diagonal elements of the first NB rows and columns of
    //          the reduced matrix.
  
    //  TAUQ    (output) DOUBLE PRECISION array dimension (NB)
    //          The scalar factors of the elementary reflectors which
    //          represent the orthogonal matrix Q. See Further Details.
  
    //  TAUP    (output) DOUBLE PRECISION array, dimension (NB)
    //          The scalar factors of the elementary reflectors which
    //          represent the orthogonal matrix P. See Further Details.
  
    //  X       (output) DOUBLE PRECISION array, dimension (LDX,NB)
    //          The m-by-nb matrix X required to update the unreduced part
    //          of A.
  
    //  LDX     (input) INTEGER
    //          The leading dimension of the array X. LDX >= M.
  
    //  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
    //          The n-by-nb matrix Y required to update the unreduced part
    //          of A.
  
    //  LDY     (output) INTEGER
    //          The leading dimension of the array Y. LDY >= N.
  
    //  Further Details
    //  ===============
  
    //  The matrices Q and P are represented as products of elementary
    //  reflectors:
  
    //     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
  
    //  Each H(i) and G(i) has the form:
  
    //     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  
    //  where tauq and taup are real scalars, and v and u are real vectors.
  
    //  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
    //  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
    //  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
  
    //  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
    //  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
    //  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
  
    //  The elements of the vectors v and u together form the m-by-nb matrix
    //  V and the nb-by-n matrix U' which are needed, with X and Y, to apply
    //  the transformation to the unreduced part of the matrix, using a block
    //  update of the form:  A := A - V*Y' - X*U'.
  
    //  The contents of A on exit are illustrated by the following examples
    //  with nb = 2:
  
    //  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
  
    //    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
    //    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
    //    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
    //    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
    //    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
    //    (  v1  v2  a   a   a  )
  
    //  where a denotes an element of the original matrix which is unchanged,
    //  vi denotes an element of the vector defining H(i), and ui an element
    //  of the vector defining G(i).
  
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
  
    //     Quick return if possible
  
    if( m <= 0 || n <= 0 ) 
        return;
  
    if( m >= n ) { 
    
        //        Reduce to upper bidiagonal form
    
        for( i = 1, i_ = i - 1, _do0 = nb; i <= _do0; i++, i_++ ) { 
      
            //           Update A(i:m,i)
      
            dgemv( 'N'/* No transpose */, m - i + 1, i - 1, -ONE, &A(0,i_), 
                   lda, &Y(0,i_), ldy, ONE, &A(i_,i_), 1 );
            dgemv( 'N'/* No transpose */, m - i + 1, i - 1, -ONE, &X(0,i_), 
                   ldx, &A(i_,0), 1, ONE, &A(i_,i_), 1 );
      
            //           Generate reflection Q(i) to annihilate A(i+1:m,i)
      
            dlarfg( m - i + 1, A(i_,i_), &A(i_,min( i + 1, m ) - 1), 
                    1, tauq[i_] );
            d[i_] = A(i_,i_);
            if( i < n ) { 
                A(i_,i_) = ONE;
        
                //              Compute Y(i+1:n,i)
        
                dgemv( 'T'/* Transpose */, m - i + 1, n - i, ONE, &A(i_ + 1,i_), 
                       lda, &A(i_,i_), 1, ZERO, &Y(i_,i_ + 1), 1 );
                dgemv( 'T'/* Transpose */, m - i + 1, i - 1, ONE, &A(0,i_), 
                       lda, &A(i_,i_), 1, ZERO, &Y(i_,0), 1 );
                dgemv( 'N'/* No transpose */, n - i, i - 1, -ONE, &Y(0,i_ + 1), 
                       ldy, &Y(i_,0), 1, ONE, &Y(i_,i_ + 1), 1 );
                dgemv( 'T'/* Transpose */, m - i + 1, i - 1, ONE, &X(0,i_), 
                       ldx, &A(i_,i_), 1, ZERO, &Y(i_,0), 1 );
                dgemv( 'T'/* Transpose */, i - 1, n - i, -ONE, &A(i_ + 1,0), 
                       lda, &Y(i_,0), 1, ONE, &Y(i_,i_ + 1), 1 );
                dscal( n - i, tauq[i_], &Y(i_,i_ + 1), 1 );
        
                //              Update A(i,i+1:n)
        
                dgemv( 'N'/* No transpose */, n - i, i, -ONE, &Y(0,i_ + 1), 
                       ldy, &A(0,i_), lda, ONE, &A(i_ + 1,i_), lda );
                dgemv( 'T'/* Transpose */, i - 1, n - i, -ONE, &A(i_ + 1,0), 
                       lda, &X(0,i_), ldx, ONE, &A(i_ + 1,i_), lda );
        
                //              Generate reflection P(i) to annihilate A(i,i+2:n)
        
                dlarfg( n - i, A(i_ + 1,i_), &A(min( i + 2, n ) - 1,i_), 
                        lda, taup[i_] );
                e[i_] = A(i_ + 1,i_);
                A(i_ + 1,i_) = ONE;
        
                //              Compute X(i+1:m,i)
        
                dgemv( 'N'/* No transpose */, m - i, n - i, ONE, &A(i_ + 1,i_ + 1), 
                       lda, &A(i_ + 1,i_), lda, ZERO, &X(i_,i_ + 1), 1 );
                dgemv( 'T'/* Transpose */, n - i, i, ONE, &Y(0,i_ + 1), 
                       ldy, &A(i_ + 1,i_), lda, ZERO, &X(i_,0), 1 );
                dgemv( 'N'/* No transpose */, m - i, i, -ONE, &A(0,i_ + 1), 
                       lda, &X(i_,0), 1, ONE, &X(i_,i_ + 1), 1 );
                dgemv( 'N'/* No transpose */, i - 1, n - i, ONE, &A(i_ + 1,0), 
                       lda, &A(i_ + 1,i_), lda, ZERO, &X(i_,0), 1 );
                dgemv( 'N'/* No transpose */, m - i, i - 1, -ONE, &X(0,i_ + 1), 
                       ldx, &X(i_,0), 1, ONE, &X(i_,i_ + 1), 1 );
                dscal( m - i, taup[i_], &X(i_,i_ + 1), 1 );
            }
        }
    }
    else { 
    
        //        Reduce to lower bidiagonal form
    
        for( i = 1, i_ = i - 1, _do1 = nb; i <= _do1; i++, i_++ ) { 
      
            //           Update A(i,i:n)
      
            dgemv( 'N'/* No transpose */, n - i + 1, i - 1, -ONE, &Y(0,i_), 
                   ldy, &A(0,i_), lda, ONE, &A(i_,i_), lda );
            dgemv( 'T'/* Transpose */, i - 1, n - i + 1, -ONE, &A(i_,0), 
                   lda, &X(0,i_), ldx, ONE, &A(i_,i_), lda );
      
            //           Generate reflection P(i) to annihilate A(i,i+1:n)
      
            dlarfg( n - i + 1, A(i_,i_), &A(min( i + 1, n ) - 1,i_), 
                    lda, taup[i_] );
            d[i_] = A(i_,i_);
            if( i < m ) { 
                A(i_,i_) = ONE;
        
                //              Compute X(i+1:m,i)
        
                dgemv( 'N'/* No transpose */, m - i, n - i + 1, ONE, 
                       &A(i_,i_ + 1), lda, &A(i_,i_), lda, ZERO, &X(i_,i_ + 1), 
                       1 );
                dgemv( 'T'/* Transpose */, n - i + 1, i - 1, ONE, &Y(0,i_), 
                       ldy, &A(i_,i_), lda, ZERO, &X(i_,0), 1 );
                dgemv( 'N'/* No transpose */, m - i, i - 1, -ONE, &A(0,i_ + 1), 
                       lda, &X(i_,0), 1, ONE, &X(i_,i_ + 1), 1 );
                dgemv( 'N'/* No transpose */, i - 1, n - i + 1, ONE, 
                       &A(i_,0), lda, &A(i_,i_), lda, ZERO, &X(i_,0), 1 );
                dgemv( 'N'/* No transpose */, m - i, i - 1, -ONE, &X(0,i_ + 1), 
                       ldx, &X(i_,0), 1, ONE, &X(i_,i_ + 1), 1 );
                dscal( m - i, taup[i_], &X(i_,i_ + 1), 1 );
        
                //              Update A(i+1:m,i)
        
                dgemv( 'N'/* No transpose */, m - i, i - 1, -ONE, &A(0,i_ + 1), 
                       lda, &Y(0,i_), ldy, ONE, &A(i_,i_ + 1), 1 );
                dgemv( 'N'/* No transpose */, m - i, i, -ONE, &X(0,i_ + 1), 
                       ldx, &A(i_,0), 1, ONE, &A(i_,i_ + 1), 1 );
        
                //              Generate reflection Q(i) to annihilate A(i+2:m,i)
        
                dlarfg( m - i, A(i_,i_ + 1), &A(i_,min( i + 2, m ) - 1), 
                        1, tauq[i_] );
                e[i_] = A(i_,i_ + 1);
                A(i_,i_ + 1) = ONE;
        
                //              Compute Y(i+1:n,i)
        
                dgemv( 'T'/* Transpose */, m - i, n - i, ONE, &A(i_ + 1,i_ + 1), 
                       lda, &A(i_,i_ + 1), 1, ZERO, &Y(i_,i_ + 1), 1 );
                dgemv( 'T'/* Transpose */, m - i, i - 1, ONE, &A(0,i_ + 1), 
                       lda, &A(i_,i_ + 1), 1, ZERO, &Y(i_,0), 1 );
                dgemv( 'N'/* No transpose */, n - i, i - 1, -ONE, &Y(0,i_ + 1), 
                       ldy, &Y(i_,0), 1, ONE, &Y(i_,i_ + 1), 1 );
                dgemv( 'T'/* Transpose */, m - i, i, ONE, &X(0,i_ + 1), 
                       ldx, &A(i_,i_ + 1), 1, ZERO, &Y(i_,0), 1 );
                dgemv( 'T'/* Transpose */, i, n - i, -ONE, &A(i_ + 1,0), 
                       lda, &Y(i_,0), 1, ONE, &Y(i_,i_ + 1), 1 );
                dscal( n - i, tauq[i_], &Y(i_,i_ + 1), 1 );
            }
        }
    }
    return;
  
    //     End of DLABRD
  
#undef  Y
#undef  X
#undef  A
} // end of function 

