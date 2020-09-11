/*
 * Default C++ implementation of dcopy
 * For optimum performance, use a machine specific bla library
 *
 * $Id: dcopy.cpp,v 1.3 1993/03/19 16:26:58 alv Exp $
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
 * $Log: dcopy.cpp,v $
 * Revision 1.3  1993/03/19 16:26:58  alv
 * added RWBLADECL linkage specification
 *
 * Revision 1.2  1993/03/05  23:07:11  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:04:34  alv
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

RWLAPKDECL void dcopy(const long &n, double dx[], const long &incx, 
                      double dy[], const long &incy)
{
    long _do0, _do1, _do2, i, i_, ix, iy, m, mp1;

  
    //     copies a vector, x, to a vector, y.
    //     uses unrolled loops for increments equal to one.
    //     jack dongarra, linpack, 3/11/78.
  
  
    if( n <= 0 ) 
        return;
    if( incx == 1 && incy == 1 ) 
        goto L_20;
  
    //        code for unequal increments or equal increments
    //          not equal to 1
  
    ix = 1;
    iy = 1;
    if( incx < 0 ) 
        ix = (-n + 1)*incx + 1;
    if( incy < 0 ) 
        iy = (-n + 1)*incy + 1;
    for( i = 1, i_ = i - 1, _do0 = n; i <= _do0; i++, i_++ ) { 
        dy[iy - 1] = dx[ix - 1];
        ix = ix + incx;
        iy = iy + incy;
    }
    return;
  
    //        code for both increments equal to 1
  
  
    //        clean-up loop
  
L_20:
    // mod() function massaged out of next line
    m = n % 7; //m = mod( n, 7 );
    if( m == 0 ) 
        goto L_40;
    for( i = 1, i_ = i - 1, _do1 = m; i <= _do1; i++, i_++ ) { 
        dy[i_] = dx[i_];
    }
    if( n < 7 ) 
        return;
L_40:
    mp1 = m + 1;
    for( i = mp1, i_ = i - 1, _do2 = n; i <= _do2; i += 7, i_ += 7 ) { 
        dy[i_] = dx[i_];
        dy[i_ + 1] = dx[i_ + 1];
        dy[i_ + 2] = dx[i_ + 2];
        dy[i_ + 3] = dx[i_ + 3];
        dy[i_ + 4] = dx[i_ + 4];
        dy[i_ + 5] = dx[i_ + 5];
        dy[i_ + 6] = dx[i_ + 6];
    }
    return;
} // end of function 

