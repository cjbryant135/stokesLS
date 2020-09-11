/*
 * C++ implementation of lapack routine dlasr
 *
 * $Id: dlasr.cpp,v 1.4 1993/04/06 20:41:27 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:12
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlasr.cpp,v $
 * Revision 1.4  1993/04/06 20:41:27  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:56  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:44  alv
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

RWLAPKDECL void /*FUNCTION*/ dlasr(const char &side, const char &pivot, const char &direct, const long &m, 
                                   const long &n, double c[], double s[], double *a, const long &lda)
{
#define A(I_,J_)  (*(a+(I_)*(lda)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, _do10, _do11, _do12, _do13, _do14, _do15, 
        _do16, _do17, _do2, _do3, _do4, _do5, _do6, _do7, _do8, _do9, 
        i, i_, info, j, j_;
    double ctemp, stemp, temp;

  
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
  
    //  DLASR   performs the transformation
  
    //     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
  
    //     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
  
    //  where A is an m by n real matrix and P is an orthogonal matrix,
    //  consisting of a sequence of plane rotations determined by the
    //  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
    //  and z = n when SIDE = 'R' or 'r' ):
  
    //  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
  
    //     P = P( z - 1 )*...*P( 2 )*P( 1 ),
  
    //  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
  
    //     P = P( 1 )*P( 2 )*...*P( z - 1 ),
  
    //  where  P( k ) is a plane rotation matrix for the following planes:
  
    //     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
    //        the plane ( k, k + 1 )
  
    //     when  PIVOT = 'T' or 't'  ( Top pivot ),
    //        the plane ( 1, k + 1 )
  
    //     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
    //        the plane ( k, z )
  
    //  c( k ) and s( k )  must contain the  cosine and sine that define the
    //  matrix  P( k ).  The two by two plane rotation part of the matrix
    //  P( k ), R( k ), is assumed to be of the form
  
    //     R( k ) = (  c( k )  s( k ) ).
    //              ( -s( k )  c( k ) )
  
    //  This version vectorises across rows of the array A when SIDE = 'L'.
  
    //  Arguments
    //  =========
  
    //  SIDE    (input) CHARACTER*1
    //          Specifies whether the plane rotation matrix P is applied to
    //          A on the left or the right.
    //          = 'L':  Left, compute A := P*A
    //          = 'R':  Right, compute A:= A*P'
  
    //  DIRECT  (input) CHARACTER*1
    //          Specifies whether P is a forward or backward sequence of
    //          plane rotations.
    //          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
    //          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
  
    //  PIVOT   (input) CHARACTER*1
    //          Specifies the plane for which P(k) is a plane rotation
    //          matrix.
    //          = 'V':  Variable pivot, the plane (k,k+1)
    //          = 'T':  Top pivot, the plane (1,k+1)
    //          = 'B':  Bottom pivot, the plane (k,z)
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix A.  If m <= 1, an immediate
    //          return is effected.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix A.  If n <= 1, an
    //          immediate return is effected.
  
    //  C, S    (input) DOUBLE PRECISION arrays, dimension
    //                  (M-1) if SIDE = 'L'
    //                  (N-1) if SIDE = 'R'
    //          c(k) and s(k) contain the cosine and sine that define the
    //          matrix P(k).  The two by two plane rotation part of the
    //          matrix P(k), R(k), is assumed to be of the form
    //          R( k ) = (  c( k )  s( k ) ).
    //                   ( -s( k )  c( k ) )
  
    //  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    //          The m by n matrix A.  On exit, A is overwritten by P*A if
    //          SIDE = 'R' or by A*P' if SIDE = 'L'.
  
    //  LDA     (input) INTEGER
    //          The leading dimension of the array A.  LDA >= max(1,M).
  
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
  
    //     Test the input parameters
  
    info = 0;
    if( !(lsame( side, 'L' ) || lsame( side, 'R' )) ) { 
        info = 1;
    }
    else if( !((lsame( pivot, 'V' ) || lsame( pivot, 'T' )) || lsame( pivot, 
                                                                      'B' )) ) { 
        info = 2;
    }
    else if( !(lsame( direct, 'F' ) || lsame( direct, 'B' )) ) { 
        info = 3;
    }
    else if( m < 0 ) { 
        info = 4;
    }
    else if( n < 0 ) { 
        info = 5;
    }
    else if( lda < max( 1, m ) ) { 
        info = 9;
    }
    if( info != 0 ) { 
        xerbla( "DLASR ", info );
        return;
    }
  
    //     Quick return if possible
  
    if( (m == 0) || (n == 0) ) 
        return;
    if( lsame( side, 'L' ) ) { 
    
        //        Form  P * A
    
        if( lsame( pivot, 'V' ) ) { 
            if( lsame( direct, 'F' ) ) { 
                for( j = 1, j_ = j - 1, _do0 = m - 1; j <= _do0; j++, j_++ ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do1 = n; i <= _do1; i++, i_++ ) { 
                            temp = A(i_,j_ + 1);
                            A(i_,j_ + 1) = ctemp*temp - stemp*A(i_,j_);
                            A(i_,j_) = stemp*temp + ctemp*A(i_,j_);
                        }
                    }
                }
            }
            else if( lsame( direct, 'B' ) ) { 
                for( j = m - 1, j_ = j - 1; j >= 1; j--, j_-- ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do2 = n; i <= _do2; i++, i_++ ) { 
                            temp = A(i_,j_ + 1);
                            A(i_,j_ + 1) = ctemp*temp - stemp*A(i_,j_);
                            A(i_,j_) = stemp*temp + ctemp*A(i_,j_);
                        }
                    }
                }
            }
        }
        else if( lsame( pivot, 'T' ) ) { 
            if( lsame( direct, 'F' ) ) { 
                for( j = 2, j_ = j - 1, _do3 = m; j <= _do3; j++, j_++ ) { 
                    ctemp = c[j_ - 1];
                    stemp = s[j_ - 1];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do4 = n; i <= _do4; i++, i_++ ) { 
                            temp = A(i_,j_);
                            A(i_,j_) = ctemp*temp - stemp*A(i_,0);
                            A(i_,0) = stemp*temp + ctemp*A(i_,0);
                        }
                    }
                }
            }
            else if( lsame( direct, 'B' ) ) { 
                for( j = m, j_ = j - 1; j >= 2; j--, j_-- ) { 
                    ctemp = c[j_ - 1];
                    stemp = s[j_ - 1];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do5 = n; i <= _do5; i++, i_++ ) { 
                            temp = A(i_,j_);
                            A(i_,j_) = ctemp*temp - stemp*A(i_,0);
                            A(i_,0) = stemp*temp + ctemp*A(i_,0);
                        }
                    }
                }
            }
        }
        else if( lsame( pivot, 'B' ) ) { 
            if( lsame( direct, 'F' ) ) { 
                for( j = 1, j_ = j - 1, _do6 = m - 1; j <= _do6; j++, j_++ ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do7 = n; i <= _do7; i++, i_++ ) { 
                            temp = A(i_,j_);
                            A(i_,j_) = stemp*A(i_,m - 1) + ctemp*temp;
                            A(i_,m - 1) = ctemp*A(i_,m - 1) - stemp*
                                temp;
                        }
                    }
                }
            }
            else if( lsame( direct, 'B' ) ) { 
                for( j = m - 1, j_ = j - 1; j >= 1; j--, j_-- ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do8 = n; i <= _do8; i++, i_++ ) { 
                            temp = A(i_,j_);
                            A(i_,j_) = stemp*A(i_,m - 1) + ctemp*temp;
                            A(i_,m - 1) = ctemp*A(i_,m - 1) - stemp*
                                temp;
                        }
                    }
                }
            }
        }
    }
    else if( lsame( side, 'R' ) ) { 
    
        //        Form A * P'
    
        if( lsame( pivot, 'V' ) ) { 
            if( lsame( direct, 'F' ) ) { 
                for( j = 1, j_ = j - 1, _do9 = n - 1; j <= _do9; j++, j_++ ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do10 = m; i <= _do10; i++, i_++ ) { 
                            temp = A(j_ + 1,i_);
                            A(j_ + 1,i_) = ctemp*temp - stemp*A(j_,i_);
                            A(j_,i_) = stemp*temp + ctemp*A(j_,i_);
                        }
                    }
                }
            }
            else if( lsame( direct, 'B' ) ) { 
                for( j = n - 1, j_ = j - 1; j >= 1; j--, j_-- ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do11 = m; i <= _do11; i++, i_++ ) { 
                            temp = A(j_ + 1,i_);
                            A(j_ + 1,i_) = ctemp*temp - stemp*A(j_,i_);
                            A(j_,i_) = stemp*temp + ctemp*A(j_,i_);
                        }
                    }
                }
            }
        }
        else if( lsame( pivot, 'T' ) ) { 
            if( lsame( direct, 'F' ) ) { 
                for( j = 2, j_ = j - 1, _do12 = n; j <= _do12; j++, j_++ ) { 
                    ctemp = c[j_ - 1];
                    stemp = s[j_ - 1];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do13 = m; i <= _do13; i++, i_++ ) { 
                            temp = A(j_,i_);
                            A(j_,i_) = ctemp*temp - stemp*A(0,i_);
                            A(0,i_) = stemp*temp + ctemp*A(0,i_);
                        }
                    }
                }
            }
            else if( lsame( direct, 'B' ) ) { 
                for( j = n, j_ = j - 1; j >= 2; j--, j_-- ) { 
                    ctemp = c[j_ - 1];
                    stemp = s[j_ - 1];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do14 = m; i <= _do14; i++, i_++ ) { 
                            temp = A(j_,i_);
                            A(j_,i_) = ctemp*temp - stemp*A(0,i_);
                            A(0,i_) = stemp*temp + ctemp*A(0,i_);
                        }
                    }
                }
            }
        }
        else if( lsame( pivot, 'B' ) ) { 
            if( lsame( direct, 'F' ) ) { 
                for( j = 1, j_ = j - 1, _do15 = n - 1; j <= _do15; j++, j_++ ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do16 = m; i <= _do16; i++, i_++ ) { 
                            temp = A(j_,i_);
                            A(j_,i_) = stemp*A(n - 1,i_) + ctemp*temp;
                            A(n - 1,i_) = ctemp*A(n - 1,i_) - stemp*
                                temp;
                        }
                    }
                }
            }
            else if( lsame( direct, 'B' ) ) { 
                for( j = n - 1, j_ = j - 1; j >= 1; j--, j_-- ) { 
                    ctemp = c[j_];
                    stemp = s[j_];
                    if( (ctemp != ONE) || (stemp != ZERO) ) { 
                        for( i = 1, i_ = i - 1, _do17 = m; i <= _do17; i++, i_++ ) { 
                            temp = A(j_,i_);
                            A(j_,i_) = stemp*A(n - 1,i_) + ctemp*temp;
                            A(n - 1,i_) = ctemp*A(n - 1,i_) - stemp*
                                temp;
                        }
                    }
                }
            }
        }
    }
  
    return;
  
    //     End of DLASR
  
#undef  A
} // end of function 

