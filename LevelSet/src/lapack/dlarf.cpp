/*
 * C++ implementation of lapack routine dlarf
 *
 * $Id: dlarf.cpp,v 1.5 1993/04/06 20:41:15 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:35:52
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlarf.cpp,v $
 * Revision 1.5  1993/04/06 20:41:15  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  18:41:23  alv
 * now passes chars explicitly, rather than indirection of a string, to shut up SUN warnings
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:15:41  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:35  alv
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

RWLAPKDECL void /*FUNCTION*/ dlarf(const char &side, const long &m, const long &n, double v[], 
                                   const long &incv, const double &tau, double *c, const long &ldc, double work[])
{
#define C(I_,J_)  (*(c+(I_)*(ldc)+(J_)))
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations


  
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
  
    //  DLARF applies a real elementary reflector H to a real m by n matrix
    //  C, from either the left or the right. H is represented in the form
  
    //        H = I - tau * v * v'
  
    //  where tau is a real scalar and v is a real vector.
  
    //  If tau = 0, then H is taken to be the unit matrix.
  
    //  Arguments
    //  =========
  
    //  SIDE    (input) CHARACTER*1
    //          = 'L': form  H * C
    //          = 'R': form  C * H
  
    //  M       (input) INTEGER
    //          The number of rows of the matrix C.
  
    //  N       (input) INTEGER
    //          The number of columns of the matrix C.
  
    //  V       (input) DOUBLE PRECISION array, dimension
    //                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
    //                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
    //          The vector v in the representation of H. V is not used if
    //          TAU = 0.
  
    //  INCV    (input) INTEGER
    //          The increment between elements of v. INCV <> 0.
  
    //  TAU     (input) DOUBLE PRECISION
    //          The value tau in the representation of H.
  
    //  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
    //          On entry, the m by n matrix C.
    //          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
    //          or C * H if SIDE = 'R'.
  
    //  LDC     (input) INTEGER
    //          The leading dimension of the array C. LDC >= max(1,M).
  
    //  WORK    (workspace) DOUBLE PRECISION array, dimension
    //                         (N) if SIDE = 'L'
    //                      or (M) if SIDE = 'R'
  
    //  =====================================================================
  
    //     .. Parameters ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
  
    if( lsame( side, 'L' ) ) { 
    
        //        Form  H * C
    
        if( tau != ZERO ) { 
      
            //           w := C' * v
      
            dgemv( 'T'/* Transpose */, m, n, ONE, c, ldc, v, incv, 
                   ZERO, work, 1 );
      
            //           C := C - v * w'
      
            dger( m, n, -tau, v, incv, work, 1, c, ldc );
        }
    }
    else { 
    
        //        Form  C * H
    
        if( tau != ZERO ) { 
      
            //           w := C * v
      
            dgemv( 'N'/* No transpose */, m, n, ONE, c, ldc, v, incv, 
                   ZERO, work, 1 );
      
            //           C := C - w * v'
      
            dger( m, n, -tau, work, 1, v, incv, c, ldc );
        }
    }
    return;
  
    //     End of DLARF
  
#undef  C
} // end of function 

