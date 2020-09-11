/*
 * C++ implementation of lapack routine dlapy2
 *
 * $Id: dlapy2.cpp,v 1.4 1993/04/06 20:41:09 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:35:40
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlapy2.cpp,v $
 * Revision 1.4  1993/04/06 20:41:09  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:33  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:25  alv
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

RWLAPKDECL double /*FUNCTION*/ dlapy2(const double &x, const double &y)
{
// PARAMETER translations
    const double ZERO = 0.0e0;
    const double ONE = 1.0e0;
// end of PARAMETER translations

    double dlapy2_v, w, xabs, yabs, z;

  
    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992
  
    //     .. Scalar Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
    //  overflow.
  
    //  Arguments
    //  =========
  
    //  X       (input) DOUBLE PRECISION
    //  Y       (input) DOUBLE PRECISION
    //          X and Y specify the values x and y.
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    xabs = fabs( x );
    yabs = fabs( y );
    w = max( xabs, yabs );
    z = min( xabs, yabs );
    if( z == ZERO ) { 
        dlapy2_v = w;
    }
    else { 
        dlapy2_v = w*sqrt( ONE + pow(z/w, 2) );
    }
    return( dlapy2_v );
  
    //     End of DLAPY2
  
} // end of function 

