/*
 * C++ implementation of lapack routine dlassq
 *
 * $Id: dlassq.cpp,v 1.6 1993/04/06 20:41:28 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:14
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlassq.cpp,v $
 * Revision 1.6  1993/04/06 20:41:28  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.5  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.4  1993/03/19  16:57:26  alv
 * sprinkled in some const
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:15:57  alv
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

RWLAPKDECL void /*FUNCTION*/ dlassq(const long &n, double x[], const long &incx, 
                                    double &scale, double &sumsq)
{
// PARAMETER translations
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, ix, ix_;
    double absxi;

  
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
  
    //  DLASSQ  returns the values  scl  and  smsq  such that
  
    //     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
  
    //  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
    //  assumed to be non-negative and  scl  returns the value
  
    //     scl = max( scale, abs( x( i ) ) ).
  
    //  scale and sumsq must be supplied in SCALE and SUMSQ and
    //  scl and smsq are overwritten on SCALE and SUMSQ respectively.
  
    //  The routine makes only one pass through the vector x.
  
    //  Arguments
    //  =========
  
    //  N       (input) INTEGER
    //          The number of elements to be used from the vector X.
  
    //  X       (input) DOUBLE PRECISION
    //          The vector for which a scaled sum of squares is computed.
    //             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
  
    //  INCX    (input) INTEGER
    //          The increment between successive values of the vector X.
    //          INCX > 0.
  
    //  SCALE   (input/output) DOUBLE PRECISION
    //          On entry, the value  scale  in the equation above.
    //          On exit, SCALE is overwritten with  scl , the scaling factor
    //          for the sum of squares.
  
    //  SUMSQ   (input/output) DOUBLE PRECISION
    //          On entry, the value  sumsq  in the equation above.
    //          On exit, SUMSQ is overwritten with  smsq , the basic sum of
    //          squares from which  scl  has been factored out.
  
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    if( n > 0 ) { 
        for( ix = 1, ix_ = ix - 1, _do0=docnt(ix,1 + (n - 1)*incx,_do1 = incx); _do0 > 0; ix += _do1, ix_ += _do1, _do0-- ) { 
            if( x[ix_] != ZERO ) { 
                absxi = fabs( x[ix_] );
                if( scale < absxi ) { 
                    sumsq = 1 + sumsq*pow(scale/absxi, 2);
                    scale = absxi;
                }
                else { 
                    sumsq = sumsq + pow(absxi/scale, 2);
                }
            }
        }
    }
    return;
  
    //     End of DLASSQ
  
} // end of function 

