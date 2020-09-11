/*
 * C++ implementation of lapack routine dorgqr
 *
 * $Id: dorgqr.cpp,v 1.6 1993/04/06 20:41:41 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:46
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dorgqr.cpp,v $
 * Revision 1.6  1993/04/06 20:41:41  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.5  1993/03/19  18:41:23  alv
 * now passes chars explicitly, rather than indirection of a string, to shut up SUN warnings
 *
 * Revision 1.4  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:16:15  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:08:00  alv
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

RWLAPKDECL void /*FUNCTION*/ dorgqr(const long &m, const long &n, const long &k, double *a, 
                                    const long &lda, double tau[], double work[], const long &lwork, long &info)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, _do2, _do3, _do4, _do5, i, i_, ib, iinfo, 
        iws, j, j_, ki, kk, l, l_, ldwork, nb, nbmin, nx;

  
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
  
    //  DORGQR generates an m by n real matrix Q with orthonormal columns,
    //  which is defined as the first n columns of a product of k elementary
    //  reflectors of order m
  
    //        Q  =  H(1) H(2) . . . H(k)
  
    //  as returned by DGEQRF.
  
    //  Arguments
    //  =========
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix Q. M >= 0.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix Q. M >= N >= 0.
  
    //  K       (input) INTEGER
    //          The number of elementary reflectors whose product defines the
    //          matrix Q. N >= K >= 0.
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          On entry, the i-th column must contain the vector which
    //          defines the elementary reflector H(i), for i = 1,2,...,k, as
    //          returned by DGEQRF in the first k columns of its array
    //          argument A.
    //          On exit, the m by n matrix Q.
  
    //  LDA     (input) INTEGER
    //          The first dimension of the array A. LDA >= max(1,M).
  
    //  TAU     (input) DOUBLE PRECISION array, dimension (K)
    //          TAU(i) must contain the scalar factor of the elementary
    //          reflector H(i), as returned by DGEQRF.
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
    //          On exit, if INFO = 0, WORK(1) returns the minimum value of
    //          LWORK required to use the optimal blocksize.
  
    //  LWORK   (input) INTEGER
    //          The dimension of the array WORK. LWORK >= max(1,N).
    //          For optimum performance LWORK should be at least N*NB, where
    //          NB is the optimal blocksize.
  
    //  INFO    (output) INTEGER
    //          = 0: successful exit
    //          < 0: if INFO = -i, the i-th argument has an illegal value
  
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
  
    //     Test the input arguments
  
    info = 0;
    if( m < 0 ) { 
        info = -1;
    }
    else if( n < 0 || n > m ) { 
        info = -2;
    }
    else if( k < 0 || k > n ) { 
        info = -3;
    }
    else if( lda < max( 1, m ) ) { 
        info = -5;
    }
    else if( lwork < max( 1, n ) ) { 
        info = -8;
    }
    if( info != 0 ) { 
        xerbla( "DORGQR", -info );
        return;
    }
  
    //     Quick return if possible
  
    if( n <= 0 ) { 
        work[0] = 1;
        return;
    }
  
    //     Determine the block size.
  
    nb = ilaenv( 1, "DORGQR", " ", m, n, k, -1 );
    nbmin = 2;
    nx = 0;
    iws = n;
    if( nb > 1 && nb < k ) { 
    
        //        Determine when to cross over from blocked to unblocked code.
    
        nx = max( 0, ilaenv( 3, "DORGQR", " ", m, n, k, -1 ) );
        if( nx < k ) { 
      
            //           Determine if workspace is large enough for blocked code.
      
            ldwork = n;
            iws = ldwork*nb;
            if( lwork < iws ) { 
        
                //              Not enough workspace to use optimal NB:  reduce NB and
                //              determine the minimum value of NB.
        
                nb = lwork/ldwork;
                nbmin = max( 2, ilaenv( 2, "DORGQR", " ", m, n, k, 
                                        -1 ) );
            }
        }
    }
  
    if( (nb >= nbmin && nb < k) && nx < k ) { 
    
        //        Use blocked code after the last block.
        //        The first kk columns are handled by the block method.
    
        ki = ((k - nx - 1)/nb)*nb;
        kk = min( k, ki + nb );
    
        //        Set A(1:kk,kk+1:n) to zero.
    
        for( j = kk + 1, j_ = j - 1, _do0 = n; j <= _do0; j++, j_++ ) { 
            for( i = 1, i_ = i - 1, _do1 = kk; i <= _do1; i++, i_++ ) { 
                A(j_,i_) = ZERO;
            }
        }
    }
    else { 
        kk = 0;
    }
  
    //     Use unblocked code for the last or only block.
  
    if( kk < n ) 
        dorg2r( m - kk, n - kk, k - kk, &A(kk,kk), lda, &tau[kk], 
                work, iinfo );
  
    if( kk > 0 ) { 
    
        //        Use blocked code
    
        for( i = ki + 1, i_ = i - 1, _do2=docnt(i,1,_do3 = -nb); _do2 > 0; i += _do3, i_ += _do3, _do2-- ) { 
            ib = min( nb, k - i + 1 );
            if( i + ib <= n ) { 
        
                //              Form the triangular factor of the block reflector
                //              H = H(i) H(i+1) . . . H(i+ib-1)
        
                dlarft( 'F'/* Forward */, 'C'/* Columnwise */, m - i + 
                        1, ib, &A(i_,i_), lda, &tau[i_], work, ldwork );
        
                //              Apply H to A(i:m,i+ib:n) from the left
        
                dlarfb( 'L'/* Left */, 'N'/* No transpose */, 'F'/* Forward */
                        , 'C'/* Columnwise */, m - i + 1, n - i - ib + 1, 
                        ib, &A(i_,i_), lda, work, ldwork, &A(i_ + ib,i_), 
                        lda, &work[ib], ldwork );
            }
      
            //           Apply H to rows i:m of current block
      
            dorg2r( m - i + 1, ib, ib, &A(i_,i_), lda, &tau[i_], work, 
                    iinfo );
      
            //           Set rows 1:i-1 of current block to zero
      
            for( j = i, j_ = j - 1, _do4 = i + ib - 1; j <= _do4; j++, j_++ ) { 
                for( l = 1, l_ = l - 1, _do5 = i - 1; l <= _do5; l++, l_++ ) { 
                    A(j_,l_) = ZERO;
                }
            }
        }
    }
  
    work[0] = iws;
    return;
  
    //     End of DORGQR
  
#undef  A
} // end of function 

