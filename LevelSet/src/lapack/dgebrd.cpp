/*
 * C++ implementation of lapack routine dgebrd
 *
 * $$RW_INSERT_HEADER "mathyrs.str"
 *
 * Translated from the Fortran using Cobalt Blue's FOR_C++,
 * and then massaged slightly to Rogue Wave format.
 *
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:33:52
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
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

RWLAPKDECL void /*FUNCTION*/ dgebrd(const long &m, const long &n, double *a, const long &lda, 
                                    double d[], double e[], double tauq[], double taup[], double work[], 
                                    const long &lwork, long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
// end of PARAMETER translations

    long _do0, _do1, _do2, _do3, i_, iinfo, j, j_, ldwrkx, 
        ldwrky, minmn, nb, nbmin, nx;
    long i=0;
    double ws;

  
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
  
    //  DGEBRD reduces a real general m by n matrix A to upper or lower
    //  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.
  
    //  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
  
    //  Arguments
    //  =========
  
    //  M       (input) INTEGER
    //          The number of rows in the matrix A.  M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns in the matrix A.  N >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the m by n general matrix to be reduced.
    //          On exit,
    //          if m >= n, the diagonal and the first superdiagonal are
    //            overwritten with the upper bidiagonal matrix B; the
    //            elements below the diagonal, with the array TAUQ, represent
    //            the orthogonal matrix Q as a product of elementary
    //            reflectors, and the elements above the first superdiagonal,
    //            with the array TAUP, represent the orthogonal matrix P as
    //            a product of elementary reflectors;
    //          if m < n, the diagonal and the first subdiagonal are
    //            overwritten with the lower bidiagonal matrix B; the
    //            elements below the first subdiagonal, with the array TAUQ,
    //            represent the orthogonal matrix Q as a product of
    //            elementary reflectors, and the elements above the diagonal,
    //            with the array TAUP, represent the orthogonal matrix P as
    //            a product of elementary reflectors.
    //          See Further Details.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
    //  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
    //          The diagonal elements of the bidiagonal matrix B:
    //          D(i) = A(i,i).
  
    //  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
    //          The off-diagonal elements of the bidiagonal matrix B:
    //          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
    //          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
  
    //  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
    //          The scalar factors of the elementary reflectors which
    //          represent the orthogonal matrix Q. See Further Details.
  
    //  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
    //          The scalar factors of the elementary reflectors which
    //          represent the orthogonal matrix P. See Further Details.
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
    //          On exit, if INFO = 0, WORK(1) returns the minimum value of
    //          LWORK required to use the optimal blocksize.
  
    //  LWORK   (input) INTEGER
    //          The length of the array WORK.  LWORK >= max(1,M,N).
    //          For optimum performance LWORK should be at least (M+N)*NB,
    //          where NB is the optimal blocksize.
  
    //  INFO    (output) INTEGER
    //          = 0: successful exit
    //          < 0: if INFO = -i, the i-th argument had an illegal value.
  
    //  Further Details
    //  ===============
  
    //  The matrices Q and P are represented as products of elementary
    //  reflectors:
  
    //  If m >= n,
  
    //     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
  
    //  Each H(i) and G(i) has the form:
  
    //     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  
    //  where tauq and taup are real scalars, and v and u are real vectors;
    //  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
    //  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
    //  tauq is stored in TAUQ(i) and taup in TAUP(i).
  
    //  If m < n,
  
    //     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
  
    //  Each H(i) and G(i) has the form:
  
    //     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  
    //  where tauq and taup are real scalars, and v and u are real vectors;
    //  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
    //  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
    //  tauq is stored in TAUQ(i) and taup in TAUP(i).
  
    //  The contents of A on exit are illustrated by the following examples:
  
    //  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
  
    //    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
    //    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
    //    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
    //    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
    //    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
    //    (  v1  v2  v3  v4  v5 )
  
    //  where d and e denote diagonal and off-diagonal elements of B, vi
    //  denotes an element of the vector defining H(i), and ui an element of
    //  the vector defining G(i).
  
    //  =====================================================================
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test the input parameters
  
    info = 0;
    if( m < 0 ) { 
        info = -1;
    }
    else if( n < 0 ) { 
        info = -2;
    }
    else if( lda < max( 1, m ) ) { 
        info = -4;
    }
    else if( lwork < vmax( 1, m, n, IEND ) ) { 
        info = -10;
    }
    if( info < 0 ) { 
        xerbla( "DGEBRD", -info );
        return;
    }
  
    //     Quick return if possible
  
    minmn = min( m, n );
    if( minmn == 0 ) { 
        work[0] = 1;
        return;
    }
  
    ws = max( m, n );
    ldwrkx = m;
    ldwrky = n;
  
    //     Set the block size NB and the crossover point NX.
  
    nb = max( 1, ilaenv( 1, "DGEBRD", " ", m, n, -1, -1 ) );
  
    if( nb > 1 && nb < minmn ) { 
    
        //        Determine when to switch from blocked to unblocked code.
    
        nx = max( nb, ilaenv( 3, "DGEBRD", " ", m, n, -1, -1 ) );
        if( nx < minmn ) { 
            ws = (m + n)*nb;
            if( lwork < ws ) { 
        
                //              Not enough work space for the optimal NB, consider using
                //              a smaller block size.
        
                nbmin = ilaenv( 2, "DGEBRD", " ", m, n, -1, -1 );
                if( lwork >= (m + n)*nbmin ) { 
                    nb = lwork/(m + n);
                }
                else { 
                    nb = 1;
                    nx = minmn;
                }
            }
        }
    }
    else { 
        nx = minmn;
    }
  
    for( i = 1, i_ = i - 1, _do0=docnt(i,minmn - nx,_do1 = nb); _do0 > 0; i += _do1, i_ += _do1, _do0-- ) { 
    
        //        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
        //        the matrices X and Y which are needed to update the unreduced
        //        part of the matrix
    
        dlabrd( m - i + 1, n - i + 1, nb, &A(i_,i_), lda, &d[i_], 
                &e[i_], &tauq[i_], &taup[i_], work, ldwrkx, &work[ldwrkx*nb], 
                ldwrky );
    
        //        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
        //        of the form  A := A - V*Y' - X*U'
    
        dgemm( 'N'/* No transpose */, 'T'/* Transpose */, m - i - nb + 
               1, n - i - nb + 1, nb, -ONE, &A(i_,i_ + nb), lda, &work[ldwrkx*nb + nb], 
               ldwrky, ONE, &A(i_ + nb,i_ + nb), lda );
        dgemm( 'N'/* No transpose */, 'N'/* No transpose */, m - i - 
               nb + 1, n - i - nb + 1, nb, -ONE, &work[nb], ldwrkx, &A(i_ + nb,i_), 
               lda, ONE, &A(i_ + nb,i_ + nb), lda );
    
        //        Copy diagonal and off-diagonal elements of B back into A
    
        if( m >= n ) { 
            for( j = i, j_ = j - 1, _do2 = i + nb - 1; j <= _do2; j++, j_++ ) { 
                A(j_,j_) = d[j_];
                A(j_ + 1,j_) = e[j_];
            }
        }
        else { 
            for( j = i, j_ = j - 1, _do3 = i + nb - 1; j <= _do3; j++, j_++ ) { 
                A(j_,j_) = d[j_];
                A(j_,j_ + 1) = e[j_];
            }
        }
    }
  
    //     Use unblocked code to reduce the remainder of the matrix
  
    dgebd2( m - i + 1, n - i + 1, &A(i - 1,i - 1), lda, &d[i - 1], 
            &e[i - 1], &tauq[i - 1], &taup[i - 1], work, iinfo );
    work[0] = ws;
    return;
  
    //     End of DGEBRD
  
#undef  A
} // end of function 

