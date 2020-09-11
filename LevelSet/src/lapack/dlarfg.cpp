/*
 * C++ implementation of lapack routine dlarfg
 *
 * $Id: dlarfg.cpp,v 1.7 1993/04/06 20:41:16 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:35:56
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlarfg.cpp,v $
 * Revision 1.7  1993/04/06 20:41:16  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.6  1993/03/19  19:35:18  alv
 * fixed constness of declaration
 *
 * Revision 1.5  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.4  1993/03/19  16:57:25  alv
 * sprinkled in some const
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:15:43  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:36  alv
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

RWLAPKDECL void /*FUNCTION*/ dlarfg(const long &n, double &alpha, double x[], const long &incx, 
                                    double &tau)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, j, j_, knt;
    double beta, rsafmn, safmin, xnorm;

  
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
  
    //  DLARFG generates a real elementary reflector H of order n, such
    //  that
  
    //        H * ( alpha ) = ( beta ),   H' * H = I.
    //            (   x   )   (   0  )
  
    //  where alpha and beta are scalars, and x is an (n-1)-element real
    //  vector. H is represented in the form
  
    //        H = I - tau * ( 1 ) * ( 1 v' ) ,
    //                      ( v )
  
    //  where tau is a real scalar and v is a real (n-1)-element
    //  vector.
  
    //  If the elements of x are all zero, then tau = 0 and H is taken to be
    //  the unit matrix.
  
    //  Otherwise  1 <= tau <= 2.
  
    //  Arguments
    //  =========
  
    //  N       (input) INTEGER
    //          The order of the elementary reflector.
  
    //  ALPHA   (input/output) DOUBLE PRECISION
    //          On entry, the value alpha.
    //          On exit, it is overwritten with the value beta.
  
    //  X       (input/output) DOUBLE PRECISION array, dimension
    //                         (1+(N-2)*abs(INCX))
    //          On entry, the vector x.
    //          On exit, it is overwritten with the vector v.
  
    //  INCX    (input) INTEGER
    //          The increment between elements of X. INCX <> 0.
  
    //  TAU     (output) DOUBLE PRECISION
    //          The value tau.
  
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
  
    if( n <= 1 ) { 
        tau = ZERO;
        return;
    }
  
    xnorm = dnrm2( n - 1, x, incx );
  
    if( xnorm == ZERO ) { 
    
        //        H  =  I
    
        tau = ZERO;
    }
    else { 
    
        //        general case
    
        beta = -sign( dlapy2( alpha, xnorm ), alpha );
        safmin = dlamch( 'S' );
        if( fabs( beta ) < safmin ) { 
      
            //           XNORM, BETA may be inaccurate; scale X and recompute them
      
            rsafmn = ONE/safmin;
            knt = 0;
        L_10:
            ;
            knt = knt + 1;
            dscal( n - 1, rsafmn, x, incx );
            beta = beta*rsafmn;
            alpha = alpha*rsafmn;
            if( fabs( beta ) < safmin ) 
                goto L_10;
      
            //           New BETA is at most 1, at least SAFMIN
      
            xnorm = dnrm2( n - 1, x, incx );
            beta = -sign( dlapy2( alpha, xnorm ), alpha );
            tau = (beta - alpha)/beta;
            dscal( n - 1, ONE/(alpha - beta), x, incx );
      
            //           If ALPHA is subnormal, it may lose relative accuracy
      
            alpha = beta;
            for( j = 1, j_ = j - 1, _do0 = knt; j <= _do0; j++, j_++ ) { 
                alpha = alpha*safmin;
            }
        }
        else { 
            tau = (beta - alpha)/beta;
            dscal( n - 1, ONE/(alpha - beta), x, incx );
            alpha = beta;
        }
    }
  
    return;
  
    //     End of DLARFG
  
} // end of function 

