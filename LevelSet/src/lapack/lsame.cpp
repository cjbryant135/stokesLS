/*
 * C++ implementation of lapack routine lsame
 *
 * $Id: lsame.cpp,v 1.4 1993/04/06 20:42:57 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:27:40
 * FOR_C++ Options SET: alloc do=rt no=p pf=xlapack s=dv str=l - prototypes
 *
 * $Log: lsame.cpp,v $
 * Revision 1.4  1993/04/06 20:42:57  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:18:03  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:09:41  alv
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

RWLAPKDECL int /*FUNCTION*/ lsame(const char &ca, const char &cb)
{
    int lsame_v;
    long inta, intb, zcode;

  
    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992
  
    //     .. Scalar Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  LSAME returns .TRUE. if CA is the same letter as CB regardless of
    //  case.
  
    //  Arguments
    //  =========
  
    //  CA      (input) CHARACTER*1
    //  CB      (input) CHARACTER*1
    //          CA and CB specify the single characters to be compared.
  
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Executable Statements ..
  
    //     Test if the characters are equal
  
    lsame_v = ca == cb;
    if( lsame_v ) 
        return( lsame_v );
  
    //     Now test for equivalence if both characters are alphabetic.
  
    zcode = 'Z';
  
    //     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
    //     machines, on which ICHAR returns a value with bit 8 set.
    //     ICHAR('A') on Prime machines returns 193 which is the same as
    //     ICHAR('A') on an EBCDIC machine.
  
    inta = ( ca );
    intb = ( cb );
  
    if( zcode == 90 || zcode == 122 ) { 
    
        //        ASCII is assumed - ZCODE is the ASCII code of either lower or
        //        upper case 'Z'.
    
        if( inta >= 97 && inta <= 122 ) 
            inta = inta - 32;
        if( intb >= 97 && intb <= 122 ) 
            intb = intb - 32;
    
    }
    else if( zcode == 233 || zcode == 169 ) { 
    
        //        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
        //        upper case 'Z'.
    
        if( ((inta >= 129 && inta <= 137) || (inta >= 145 && inta <= 
                                              153)) || (inta >= 162 && inta <= 169) ) 
            inta = inta + 64;
        if( ((intb >= 129 && intb <= 137) || (intb >= 145 && intb <= 
                                              153)) || (intb >= 162 && intb <= 169) ) 
            intb = intb + 64;
    
    }
    else if( zcode == 218 || zcode == 250 ) { 
    
        //        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
        //        plus 128 of either lower or upper case 'Z'.
    
        if( inta >= 225 && inta <= 250 ) 
            inta = inta - 32;
        if( intb >= 225 && intb <= 250 ) 
            intb = intb - 32;
    }
    lsame_v = inta == intb;
  
    //     RETURN
  
    //     End of LSAME
  
    return( lsame_v );
} // end of function 

