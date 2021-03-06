/*
 * Default C++ implementation of xerbla
 * For optimum performance, use a machine specific bla library
 *
 * $Id: xerbla.cpp,v 1.7 1993/06/29 22:20:04 alv Exp $
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
 * Rewritten from scratch in C++ to avoid using Fortran I/O.
 * This needs to be set up to throw an exception.
 *
 * $Log: xerbla.cpp,v $
 * Revision 1.7  1993/06/29 22:20:04  alv
 * added include of lapack.h to get defn of RWLAPKDECL for math.h++ <= v5.0.0
 *
 * Revision 1.6  1993/04/06  20:42:59  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.5  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.4  1993/03/19  16:57:29  alv
 * sprinkled in some const
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:18:05  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:09:48  alv
 * Initial revision
 *
 */

#include <stdlib.h>
#include <iostream.h>
#if 0
#include "rw/lapkdefs.h"
#include "rw/bla.h"
#include "rw/fortran.h" /* Fortran run time library */
#include "rw/lapack.h"
#else
#include "level/lapack.h"
#endif

RWLAPKDECL void /*FUNCTION*/ xerbla(char *srname, const long &info)
{
    //  -- LAPACK auxiliary routine (preliminary version) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992
  
    //     .. Scalar Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  XERBLA  is an error handler for the LAPACK routines.
    //  It is called by an LAPACK routine if an input parameter has an
    //  invalid value.  A message is printed and execution stops.
  
    //  Installers may consider modifying the STOP statement in order to
    //  call system-specific exception-handling facilities.
  
    //  Arguments
    //  =========
  
    //  SRNAME  (input) CHARACTER*6
    //          The name of the routine which called XERBLA.
  
    //  INFO    (input) INTEGER
    //          The position of the invalid parameter in the parameter list
    //          of the calling routine.
  
    cerr <<    " ***\n";
    cerr <<    " *** xerbla called in routine '" << srname;
    cerr << "'\n *** error number is " << info;
    cerr <<  "\n *** aborting\n";
    cerr <<    " *** " << endl;

    exit((int)info);

} // end of function 

