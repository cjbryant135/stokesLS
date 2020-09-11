/*
 * C++ implementation of lapack routine dlamch
 *
 * $Id: dlamch.cpp,v 1.7 1993/04/06 20:41:01 alv Exp $
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
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:35:21
 * FOR_C++ Options SET: alloc do=rt no=p pf=dlapack,xlapack,dbla s=dv str=l - prototypes
 *
 * $Log: dlamch.cpp,v $
 * Revision 1.7  1993/04/06 20:41:01  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.6  1993/03/19  17:25:30  alv
 * made subroutine dlamch1,... have RWLAPKDECL linkage
 *
 * Revision 1.5  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.4  1993/03/19  16:57:23  alv
 * sprinkled in some const
 *
 * Revision 1.3  1993/03/09  16:14:40  alv
 * made parms const
 *
 * Revision 1.2  1993/03/05  23:15:23  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:07:15  alv
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

RWLAPKDECL void dlamc1(long&,long&,int&,int&);
RWLAPKDECL void dlamc2(long&,long&,int&,double&,long&,double&,long&,double&);
RWLAPKDECL double dlamc3(const double&,const double&);
RWLAPKDECL void dlamc4(long&,const double&,const long&);
RWLAPKDECL void dlamc5(const long&,const long&,const long&,const int&,long&,double&);

RWLAPKDECL double /*FUNCTION*/ dlamch(const char &cmach)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    int lrnd;
    long beta, imax, imin, it;
    double dlamch_v, rmach, small;
    static double base, emax, emin, eps, prec, rmax, rmin, rnd, sfmin,
        t;
    static int first = TRUE;

    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992

    //     .. Scalar Arguments ..
    //     ..

    //  Purpose
    //  =======

    //  DLAMCH determines double precision machine parameters.

    //  Arguments
    //  =========

    //  CMACH   (input) CHARACTER*1
    //          Specifies the value to be returned by DLAMCH:
    //          = 'E' or 'e',   DLAMCH := eps
    //          = 'S' or 's ,   DLAMCH := sfmin
    //          = 'B' or 'b',   DLAMCH := base
    //          = 'P' or 'p',   DLAMCH := eps*base
    //          = 'N' or 'n',   DLAMCH := t
    //          = 'R' or 'r',   DLAMCH := rnd
    //          = 'M' or 'm',   DLAMCH := emin
    //          = 'U' or 'u',   DLAMCH := rmin
    //          = 'L' or 'l',   DLAMCH := emax
    //          = 'O' or 'o',   DLAMCH := rmax

    //          where

    //          eps   = relative machine precision
    //          sfmin = safe minimum, such that 1/sfmin does not overflow
    //          base  = base of the machine
    //          prec  = eps*base
    //          t     = number of (base) digits in the mantissa
    //          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
    //          emin  = minimum exponent before (gradual) underflow
    //          rmin  = underflow threshold - base**(emin-1)
    //          emax  = largest exponent before overflow
    //          rmax  = overflow threshold  - (base**emax)*(1-eps)


    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..

    if( first ) {
        first = FALSE;
        dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax );
        base = beta;
        t = it;
        if( lrnd ) {
            rnd = ONE;
            eps = (pow(base, 1 - it))/2;
        }
        else {
            rnd = ZERO;
            eps = pow(base, 1 - it);
        }
        prec = eps*base;
        emin = imin;
        emax = imax;
        sfmin = rmin;
        small = ONE/rmax;
        if( small >= sfmin ) {

            //           Use SMALL plus a bit, to avoid the possibility of rounding
            //           causing overflow when computing  1/sfmin.

            sfmin = small*(ONE + eps);
        }
    }

    if( lsame( cmach, 'E' ) ) {
        rmach = eps;
    }
    else if( lsame( cmach, 'S' ) ) {
        rmach = sfmin;
    }
    else if( lsame( cmach, 'B' ) ) {
        rmach = base;
    }
    else if( lsame( cmach, 'P' ) ) {
        rmach = prec;
    }
    else if( lsame( cmach, 'N' ) ) {
        rmach = t;
    }
    else if( lsame( cmach, 'R' ) ) {
        rmach = rnd;
    }
    else if( lsame( cmach, 'M' ) ) {
        rmach = emin;
    }
    else if( lsame( cmach, 'U' ) ) {
        rmach = rmin;
    }
    else if( lsame( cmach, 'L' ) ) {
        rmach = emax;
    }
    else if( lsame( cmach, 'O' ) ) {
        rmach = rmax;
    }
  
    dlamch_v = rmach;
    return( dlamch_v );

    //     End of DLAMCH

} // end of function

// ***********************************************************************

RWLAPKDECL void /*FUNCTION*/ dlamc1(long &beta, long &t, int &rnd,
                                    int &ieee1)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    static int lieee1, lrnd;
    static long lbeta, lt;
    double a, b, c, f, one, qtr, savec, t1, t2;
    static int first = TRUE;


    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992

    //     .. Scalar Arguments ..
    //     ..

    //  Purpose
    //  =======

    //  DLAMC1 determines the machine parameters given by BETA, T, RND, and
    //  IEEE1.

    //  Arguments
    //  =========

    //  BETA    (output) INTEGER
    //          The base of the machine.

    //  T       (output) INTEGER
    //          The number of ( BETA ) digits in the mantissa.

    //  RND     (output) LOGICAL
    //          Specifies whether proper rounding  ( RND = .TRUE. )  or
    //          chopping  ( RND = .FALSE. )  occurs in addition. This may not
    //          be a reliable guide to the way in which the machine performs
    //          its arithmetic.

    //  IEEE1   (output) LOGICAL
    //          Specifies whether rounding appears to be done in the IEEE
    //          'round to nearest' style.

    //  Further Details
    //  ===============

    //  The routine is based on the routine  ENVRON  by Malcolm and
    //  incorporates suggestions by Gentleman and Marovich. See

    //     Malcolm M. A. (1972) Algorithms to reveal properties of
    //        floating-point arithmetic. Comms. of the ACM, 15, 949-951.

    //     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
    //        that reveal properties of floating point arithmetic units.
    //        Comms. of the ACM, 17, 276-277.


    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..

    if( first ) {
        first = FALSE;
        one = 1;

        //        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
        //        IEEE1, T and RND.

        //        Throughout this routine  we use the function  DLAMC3  to ensure
        //        that relevant values are  stored and not held in registers,  or
        //        are not affected by optimizers.

        //        Compute  a = 2.0**m  with the  smallest positive integer m such
        //        that

        //           fl( a + 1.0 ) = a.

        a = 1;
        c = 1;

        //+       WHILE( C.EQ.ONE )LOOP
    L_10:
        ;
        if( c == one ) {
            a = 2*a;
            c = dlamc3( a, one );
            c = dlamc3( c, -a );
            goto L_10;
        }
        //+       END WHILE

        //        Now compute  b = 2.0**m  with the smallest positive integer m
        //        such that

        //           fl( a + b ) .gt. a.

        b = 1;
        c = dlamc3( a, b );

        //+       WHILE( C.EQ.A )LOOP
    L_20:
        ;
        if( c == a ) {
            b = 2*b;
            c = dlamc3( a, b );
            goto L_20;
        }
        //+       END WHILE

        //        Now compute the base.  a and c  are neighbouring floating point
        //        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
        //        their difference is beta. Adding 0.25 to c is to ensure that it
        //        is truncated to beta and not ( beta - 1 ).

        qtr = one/4;
        savec = c;
        c = dlamc3( c, -a );
        lbeta = c + qtr;

        //        Now determine whether rounding or chopping occurs,  by adding a
        //        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.

        b = lbeta;
        f = dlamc3( b/2, -b/100 );
        c = dlamc3( f, a );
        if( c == a ) {
            lrnd = TRUE;
        }
        else {
            lrnd = FALSE;
        }
        f = dlamc3( b/2, b/100 );
        c = dlamc3( f, a );
        if( (lrnd) && (c == a) )
            lrnd = FALSE;

        //        Try and decide whether rounding is done in the  IEEE  'round to
        //        nearest' style. B/2 is half a unit in the last place of the two
        //        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
        //        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
        //        A, but adding B/2 to SAVEC should change SAVEC.

        t1 = dlamc3( b/2, a );
        t2 = dlamc3( b/2, savec );
        lieee1 = ((t1 == a) && (t2 > savec)) && lrnd;

        //        Now find  the  mantissa, t.  It should  be the  integer part of
        //        log to the base beta of a,  however it is safer to determine  t
        //        by powering.  So we find t as the smallest positive integer for
        //        which

        //           fl( beta**t + 1.0 ) = 1.0.

        lt = 0;
        a = 1;
        c = 1;

        //+       WHILE( C.EQ.ONE )LOOP
    L_30:
        ;
        if( c == one ) {
            lt = lt + 1;
            a = a*lbeta;
            c = dlamc3( a, one );
            c = dlamc3( c, -a );
            goto L_30;
        }
        //+       END WHILE

    }

    beta = lbeta;
    t = lt;
    rnd = lrnd;
    ieee1 = lieee1;
    return;

    //     End of DLAMC1

} // end of function

// ***********************************************************************

RWLAPKDECL void /*FUNCTION*/ dlamc2(long &beta, long &t, int &rnd,
                                    double &eps, long &emin, double &rmin, long &emax, double &rmax)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    int ieee, lieee1, lrnd;
    long _do0, gnmin, gpmin, i, i_, ngnmin, ngpmin;
    static long lbeta, lemax, lemin, lt;
    double a, b, c, half, one, rbase, sixth, small, third, two, zero;
    static double leps, lrmax, lrmin;
    static int first = TRUE;
    static int iwarn = FALSE;
#if 0
    static F77FMT _fmts[] = {
        9999, "(//' WARNING. The value EMIN may be incorrect:-','  EMIN = ',\
i8,/' If, after inspection, the value EMIN looks',\
' acceptable please comment out ',/\
' the IF block as marked within the code of routine',' DLAMC2,',/\
' otherwise supply EMIN explicitly.',/)",
        0L,"" };
#endif

    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992

    //     .. Scalar Arguments ..
    //     ..

    //  Purpose
    //  =======

    //  DLAMC2 determines the machine parameters specified in its argument
    //  list.

    //  Arguments
    //  =========

    //  BETA    (output) INTEGER
    //          The base of the machine.

    //  T       (output) INTEGER
    //          The number of ( BETA ) digits in the mantissa.

    //  RND     (output) LOGICAL
    //          Specifies whether proper rounding  ( RND = .TRUE. )  or
    //          chopping  ( RND = .FALSE. )  occurs in addition. This may not
    //          be a reliable guide to the way in which the machine performs
    //          its arithmetic.

    //  EPS     (output) DOUBLE PRECISION
    //          The smallest positive number such that

    //             fl( 1.0 - EPS ) .LT. 1.0,

    //          where fl denotes the computed value.

    //  EMIN    (output) INTEGER
    //          The minimum exponent before (gradual) underflow occurs.

    //  RMIN    (output) DOUBLE PRECISION
    //          The smallest normalized number for the machine, given by
    //          BASE**( EMIN - 1 ), where  BASE  is the floating point value
    //          of BETA.

    //  EMAX    (output) INTEGER
    //          The maximum exponent before overflow occurs.

    //  RMAX    (output) DOUBLE PRECISION
    //          The largest positive number for the machine, given by
    //          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
    //          value of BETA.

    //  Further Details
    //  ===============

    //  The computation of  EPS  is based on a routine PARANOIA by
    //  W. Kahan of the University of California at Berkeley.


    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..

    if( first ) {
        first = FALSE;
        zero = 0;
        one = 1;
        two = 2;

        //        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
        //        BETA, T, RND, EPS, EMIN and RMIN.

        //        Throughout this routine  we use the function  DLAMC3  to ensure
        //        that relevant values are stored  and not held in registers,  or
        //        are not affected by optimizers.

        //        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.

        dlamc1( lbeta, lt, lrnd, lieee1 );

        //        Start to find EPS.

        b = lbeta;
        a = pow(b, -lt);
        leps = a;

        //        Try some tricks to see whether or not this is the correct  EPS.

        b = two/3;
        half = one/2;
        sixth = dlamc3( b, -half );
        third = dlamc3( sixth, sixth );
        b = dlamc3( third, -half );
        b = dlamc3( b, sixth );
        b = fabs( b );
        if( b < leps )
            b = leps;

        leps = 1;

        //+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
    L_10:
        ;
        if( (leps > b) && (b > zero) ) {
            leps = b;
            c = dlamc3( half*leps, (pow(two, 5))*(pow(leps, 2)) );
            c = dlamc3( half, -c );
            b = dlamc3( half, c );
            c = dlamc3( half, -b );
            b = dlamc3( half, c );
            goto L_10;
        }
        //+       END WHILE

        if( a < leps )
            leps = a;

        //        Computation of EPS complete.

        //        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
        //        Keep dividing  A by BETA until (gradual) underflow occurs. This
        //        is detected when we cannot recover the previous A.

        rbase = one/lbeta;
        small = one;
        for( i = 1, i_ = i - 1; i <= 3; i++, i_++ ) {
            small = dlamc3( small*rbase, zero );
        }
        a = dlamc3( one, small );
        dlamc4( ngpmin, one, lbeta );
        dlamc4( ngnmin, -one, lbeta );
        dlamc4( gpmin, a, lbeta );
        dlamc4( gnmin, -a, lbeta );
        ieee = FALSE;

        if( (ngpmin == ngnmin) && (gpmin == gnmin) ) {
            if( ngpmin == gpmin ) {
                lemin = ngpmin;
                //            ( Non twos-complement machines, no gradual underflow;
                //              e.g.,  VAX )
            }
            else if( (gpmin - ngpmin) == 3 ) {
                lemin = ngpmin - 1 + lt;
                ieee = TRUE;
                //            ( Non twos-complement machines, with gradual underflow;
                //              e.g., IEEE standard followers )
            }
            else {
                lemin = min( ngpmin, gpmin );
                //            ( A guess; no known machine )
                iwarn = TRUE;
            }

        }
        else if( (ngpmin == gpmin) && (ngnmin == gnmin) ) {
            if( fabs( ngpmin - ngnmin ) == 1 ) {
                lemin = max( ngpmin, ngnmin );
                //            ( Twos-complement machines, no gradual underflow;
                //              e.g., CYBER 205 )
            }
            else {
                lemin = min( ngpmin, ngnmin );
                //            ( A guess; no known machine )
                iwarn = TRUE;
            }

        }
        else if( (fabs( ngpmin - ngnmin ) == 1) && (gpmin == gnmin) ) {
            if( (gpmin - min( ngpmin, ngnmin )) == 3 ) {
                lemin = max( ngpmin, ngnmin ) - 1 + lt;
                //            ( Twos-complement machines with gradual underflow;
                //              no known machine )
            }
            else {
                lemin = min( ngpmin, ngnmin );
                //            ( A guess; no known machine )
                iwarn = TRUE;
            }

        }
        else {
            lemin = vmin( ngpmin, ngnmin, gpmin, gnmin, IEND );
            //         ( A guess; no known machine )
            iwarn = TRUE;
        }
        // **
        // Comment out this if block if EMIN is ok
#if 0
        if( iwarn ) {
            first = TRUE;
            writef( 6, FMTR(9999), "%ld\n", lemin );
        }
#endif
        // **

        //        Assume IEEE arithmetic if we found denormalised  numbers above,
        //        or if arithmetic seems to round in the  IEEE style,  determined
        //        in routine DLAMC1. A true IEEE machine should have both  things
        //        true; however, faulty machines may have one or the other.

        ieee = ieee || lieee1;

        //        Compute  RMIN by successive division by  BETA. We could compute
        //        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
        //        this computation.

        lrmin = 1;
        for( i = 1, i_ = i - 1, _do0 = 1 - lemin; i <= _do0; i++, i_++ ) {
            lrmin = dlamc3( lrmin*rbase, zero );
        }

        //        Finally, call DLAMC5 to compute EMAX and RMAX.

        dlamc5( lbeta, lt, lemin, ieee, lemax, lrmax );
    }

    beta = lbeta;
    t = lt;
    rnd = lrnd;
    eps = leps;
    emin = lemin;
    rmin = lrmin;
    emax = lemax;
    rmax = lrmax;

    return;


    //     End of DLAMC2

} // end of function

// ***********************************************************************

RWLAPKDECL double /*FUNCTION*/ dlamc3(const double &a, const double &b)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    double dlamc3_v;


    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992

    //     .. Scalar Arguments ..
    //     ..

    //  Purpose
    //  =======

    //  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
    //  the addition of  A  and  B ,  for use in situations where optimizers
    //  might hold one of these in a register.

    //  Arguments
    //  =========

    //  A, B    (input) DOUBLE PRECISION
    //          The values A and B.


    //     .. Executable Statements ..

    dlamc3_v = a + b;

    return( dlamc3_v );

    //     End of DLAMC3

} // end of function

// ***********************************************************************

RWLAPKDECL void /*FUNCTION*/ dlamc4(long &emin,const double &start, const long &base)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

// Visual Age compiler throws fatal exception for underflows and overflows.
// This circumvents the forced underflow/overflow by hardcoding the value of emin
#ifdef __IBMCPP__
    emin = -1020;
#else
    long _do0, _do1, i, i_;
    double a, b1, b2, c1, c2, d1, d2, one, rbase, zero;


    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992

    //     .. Scalar Arguments ..
    //     ..

    //  Purpose
    //  =======

    //  DLAMC4 is a service routine for DLAMC2.

    //  Arguments
    //  =========

    //  EMIN    (output) EMIN
    //          The minimum exponent before (gradual) underflow, computed by
    //          setting A = START and dividing by BASE until the previous A
    //          can not be recovered.

    //  START   (input) DOUBLE PRECISION
    //          The starting point for determining EMIN.

    //  BASE    (input) INTEGER
    //          The base of the machine.


    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..

    a = start;
    one = 1;
    rbase = one/base;
    zero = 0;
    emin = 1;
    b1 = dlamc3( a*rbase, zero );
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
    //+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
    //    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
L_10:
    ;
    if( (((c1 == a) && (c2 == a)) && (d1 == a)) && (d2 == a) ) {
        emin = emin - 1;
        a = b1;
        b1 = dlamc3( a/base, zero );
        c1 = dlamc3( b1*base, zero );
        d1 = zero;
        for( i = 1, i_ = i - 1, _do0 = base; i <= _do0; i++, i_++ ) {
            d1 = d1 + b1;
        }
        b2 = dlamc3( a*rbase, zero );
        c2 = dlamc3( b2/rbase, zero );
        d2 = zero;
        for( i = 1, i_ = i - 1, _do1 = base; i <= _do1; i++, i_++ ) {
            d2 = d2 + b2;
        }
        goto L_10;
    }
    //+    END WHILE
#endif
    return;

    //     End of DLAMC4

} // end of function

// ***********************************************************************

RWLAPKDECL void /*FUNCTION*/ dlamc5(const long &beta, const long &p, const long &emin,
                                    const int &ieee, long &emax, double &rmax)
{
// PARAMETER translations
    const double ONE = 1.0e0;
    const double ZERO = 0.0e0;
// end of PARAMETER translations

    long _do0, _do1, exbits, expsum, i, i_, lexp, nbits, try_,
        uexp;
    double oldy, recbas, y, z;


    //  -- LAPACK auxiliary routine (version 1.0) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 29, 1992

    //     .. Scalar Arguments ..
    //     ..

    //  Purpose
    //  =======

    //  DLAMC5 attempts to compute RMAX, the largest machine floating-point
    //  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
    //  approximately to a power of 2.  It will fail on machines where this
    //  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
    //  EMAX = 28718).  It will also fail if the value supplied for EMIN is
    //  too large (i.e. too close to zero), probably with overflow.

    //  Arguments
    //  =========

    //  BETA    (input) INTEGER
    //          The base of floating-point arithmetic.

    //  P       (input) INTEGER
    //          The number of base BETA digits in the mantissa of a
    //          floating-point value.

    //  EMIN    (input) INTEGER
    //          The minimum exponent before (gradual) underflow.

    //  IEEE    (input) LOGICAL
    //          A logical flag specifying whether or not the arithmetic
    //          system is thought to comply with the IEEE standard.

    //  EMAX    (output) INTEGER
    //          The largest exponent before overflow

    //  RMAX    (output) DOUBLE PRECISION
    //          The largest machine floating-point number.


    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..

    //     First compute LEXP and UEXP, two powers of 2 that bound
    //     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
    //     approximately to the bound that is closest to abs(EMIN).
    //     (EMAX is the exponent of the required number RMAX).

    lexp = 1;
    exbits = 1;
L_10:
    ;
    try_ = lexp*2;
    if( try_ <= (-emin) ) {
        lexp = try_;
        exbits = exbits + 1;
        goto L_10;
    }
    if( lexp == -emin ) {
        uexp = lexp;
    }
    else {
        uexp = try_;
        exbits = exbits + 1;
    }

    //     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
    //     than or equal to EMIN. EXBITS is the number of bits needed to
    //     store the exponent.

    if( (uexp + emin) > (-lexp - emin) ) {
        expsum = 2*lexp;
    }
    else {
        expsum = 2*uexp;
    }

    //     EXPSUM is the exponent range, approximately equal to
    //     EMAX - EMIN + 1 .

    emax = expsum + emin - 1;
    nbits = 1 + exbits + p;

    //     NBITS is the total number of bits needed to store a
    //     floating-point number.

    if( (nbits%2/*mod( nbits, 2 )*/ == 1) && (beta == 2) ) {

        //        Either there are an odd number of bits used to store a
        //        floating-point number, which is unlikely, or some bits are
        //        not used in the representation of numbers, which is possible,
        //        (e.g. Cray machines) or the mantissa has an implicit bit,
        //        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
        //        most likely. We have to assume the last alternative.
        //        If this is true, then we need to reduce EMAX by one because
        //        there must be some way of representing zero in an implicit-bit
        //        system. On machines like Cray, we are reducing EMAX by one
        //        unnecessarily.

        emax = emax - 1;
    }

    if( ieee ) {

        //        Assume we are on an IEEE machine which reserves one exponent
        //        for infinity and NaN.

        emax = emax - 1;
    }
// Visual Age compiler throws fatal exception for underflows and overflows.
// This circumvents the forced underflow/overflow by hardcoding the value of emax
#ifdef __IBMCPP__
    emax = 1020;
#endif
    //     Now create RMAX, the largest machine number, which should
    //     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .

    //     First compute 1.0 - BETA**(-P), being careful that the
    //     result is less than 1.0 .

    recbas = ONE/beta;
    z = beta - ONE;
    y = ZERO;
    for( i = 1, i_ = i - 1, _do0 = p; i <= _do0; i++, i_++ ) {
        z = z*recbas;
        if( y < ONE )
            oldy = y;
        y = dlamc3( y, z );
    }
    if( y >= ONE )
        y = oldy;

    //     Now multiply by BETA**EMAX to get RMAX.

    for( i = 1, i_ = i - 1, _do1 = emax; i <= _do1; i++, i_++ ) {
        y = dlamc3( y*beta, ZERO );
    }

    rmax = y;
    return;

    //     End of DLAMC5

} // end of function

