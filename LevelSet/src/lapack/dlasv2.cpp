/*
 * C++ implementation of lapack routine dlasv2
 *
 * $Id: dlasv2.cpp,v 1.6 1993/04/06 20:41:28 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:36:15
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlasv2.cpp,v $
 * Revision 1.6  1993/04/06 20:41:28  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.5  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.4  1993/03/19  16:57:27  alv
 * sprinkled in some const
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:15:57  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:45  alv
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

RWLAPKDECL void /*FUNCTION*/ dlasv2(const double &f, const double &g, const double &h, double &ssmin, 
                                    double &ssmax, double &snr, double &csr, double &snl, double &csl)
{
// PARAMETER translations
    const double ZERO = 0.0e0;
    const double HALF = 0.5e0;
    const double ONE = 1.0e0;
    const double TWO = 2.0e0;
    const double FOUR = 4.0e0;
// end of PARAMETER translations

    int gasmal, swap;
    long pmax;
    double a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m, mm, r, s, 
        slt, srt, t, temp, tsign, tt;

  
    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992
  
    //     .. Scalar Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  DLASV2 computes the singular value decomposition of a 2-by-2
    //  triangular matrix
    //     [  F   G  ]
    //     [  0   H  ].
    //  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
    //  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
    //  right singular vectors for abs(SSMAX), giving the decomposition
  
    //     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
    //     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
  
    //  Arguments
    //  =========
  
    //  F       (input) DOUBLE PRECISION
    //          The (1,1) entry of the 2-by-2 matrix.
  
    //  G       (input) DOUBLE PRECISION
    //          The (1,2) entry of the 2-by-2 matrix.
  
    //  H       (input) DOUBLE PRECISION
    //          The (2,2) entry of the 2-by-2 matrix.
  
    //  SSMIN   (output) DOUBLE PRECISION
    //          abs(SSMIN) is the smaller singular value.
  
    //  SSMAX   (output) DOUBLE PRECISION
    //          abs(SSMAX) is the larger singular value.
  
    //  SNL     (output) DOUBLE PRECISION
    //  CSL     (output) DOUBLE PRECISION
    //          The vector (CSL, SNL) is a unit left singular vector for the
    //          singular value abs(SSMAX).
  
    //  SNR     (output) DOUBLE PRECISION
    //  CSR     (output) DOUBLE PRECISION
    //          The vector (CSR, SNR) is a unit right singular vector for the
    //          singular value abs(SSMAX).
  
    //  Further Details
    //  ===============
  
    //  Any input parameter may be aliased with any output parameter.
  
    //  Barring over/underflow and assuming a guard digit in subtraction, all
    //  output quantities are correct to within a few units in the last
    //  place (ulps).
  
    //  In IEEE arithmetic, the code works correctly if one matrix entry is
    //  infinite.
  
    //  Overflow will not occur unless the largest singular value itself
    //  overflows or is within a few ulps of overflow. (On machines with
    //  partial overflow, like the Cray, overflow may occur if the largest
    //  singular value is within a factor of 2 of overflow.)
  
    //  Underflow is harmless if underflow is gradual. Otherwise, results
    //  may correspond to a matrix modified by perturbations of size near
    //  the underflow threshold.
  
  
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
  
    ft = f;
    fa = fabs( ft );
    ht = h;
    ha = fabs( h );
  
    //     PMAX points to the maximum absolute entry of matrix
    //       PMAX = 1 if F largest in absolute values
    //       PMAX = 2 if G largest in absolute values
    //       PMAX = 3 if H largest in absolute values
  
    pmax = 1;
    swap = ha > fa;
    if( swap ) { 
        pmax = 3;
        temp = ft;
        ft = ht;
        ht = temp;
        temp = fa;
        fa = ha;
        ha = temp;
    
        //        Now FA .ge. HA
    
    }
    gt = g;
    ga = fabs( gt );
    if( ga == ZERO ) { 
    
        //        Diagonal matrix
    
        ssmin = ha;
        ssmax = fa;
        clt = ONE;
        crt = ONE;
        slt = ZERO;
        srt = ZERO;
    }
    else { 
        gasmal = TRUE;
        if( ga > fa ) { 
            pmax = 2;
            if( ONE + (fa/ga) == ONE ) { 
        
                //              Case of very large GA
        
                gasmal = FALSE;
                ssmax = ga;
                if( ha > ONE ) { 
                    ssmin = fa/(ga/ha);
                }
                else { 
                    ssmin = (fa/ga)*ha;
                }
                clt = ONE;
                slt = ht/gt;
                srt = ONE;
                crt = ft/gt;
            }
        }
        if( gasmal ) { 
      
            //           Normal case
      
            d = fa - ha;
            if( d == fa ) { 
        
                //              Copes with infinite F or H
        
                l = ONE;
            }
            else { 
                l = d/fa;
            }
      
            //           Note that 0 .le. L .le. 1
      
            m = gt/ft;
      
            //           Note that abs(M) .le. 1/macheps
      
            t = TWO - l;
      
            //           Note that T .ge. 1
      
            mm = m*m;
            tt = t*t;
            s = sqrt( tt + mm );
      
            //           Note that 1 .le. S .le. 1 + 1/macheps
      
            if( l == ZERO ) { 
                r = fabs( m );
            }
            else { 
                r = sqrt( l*l + mm );
            }
      
            //           Note that 0 .le. R .le. 1 + 1/macheps
      
            a = HALF*(s + r);
      
            //           Note that 1 .le. A .le. 1 + abs(M)
      
            ssmin = ha/a;
            ssmax = fa*a;
            if( mm == ZERO ) { 
        
                //              Note that M is very tiny
        
                if( l == ZERO ) { 
                    t = sign( TWO, ft )*sign( ONE, gt );
                }
                else { 
                    t = gt/sign( d, ft ) + m/t;
                }
            }
            else { 
                t = (m/(s + t) + m/(r + l))*(ONE + a);
            }
            l = sqrt( t*t + FOUR );
            crt = TWO/l;
            srt = t/l;
            clt = (crt + srt*m)/a;
            slt = (ht/ft)*srt/a;
        }
    }
    if( swap ) { 
        csl = srt;
        snl = crt;
        csr = slt;
        snr = clt;
    }
    else { 
        csl = clt;
        snl = slt;
        csr = crt;
        snr = srt;
    }
  
    //     Correct signs of SSMAX and SSMIN
  
    if( pmax == 1 ) 
        tsign = sign( ONE, csr )*sign( ONE, csl )*sign( ONE, f );
    if( pmax == 2 ) 
        tsign = sign( ONE, snr )*sign( ONE, csl )*sign( ONE, g );
    if( pmax == 3 ) 
        tsign = sign( ONE, snr )*sign( ONE, snl )*sign( ONE, h );
    ssmax = sign( ssmax, tsign );
    ssmin = sign( ssmin, tsign*sign( ONE, f )*sign( ONE, h ) );
    return;
  
    //     End of DLASV2
  
} // end of function 

