/*
 * Fortran run time library functions needed by lapack.h++
 *
 * $Id: fortran.cpp,v 1.2 1993/07/21 19:01:12 alv Exp $
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
 * $Log: fortran.cpp,v $
 * Revision 1.2  1993/07/21 19:01:12  alv
 * ported to Microsoft C++ v8
 *
 * Revision 1.1  1993/03/03  16:09:35  alv
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include <iostream>
#if 0
#include "rw/fortran.h"
#else
#include "lapack.h"
#endif

/*
 * memerr() is called if an attempt to allocate a large array  off the
 * heap fails.  This should call RWThrow.
 */

void memerr(char *s)
{
    std::cerr << "Memory allocation error in Fortran translated code\n";
    std::cerr << "Message: '" << s << "'" << std::endl;
    exit(1);
}

/*
 * pow functions
 *
 * From Cobalt Blue run time library, which used these comments:
 *  Very slightly modified version of power() from Computer Language, Sept. 86,
 *  pg 87, by Jon Snader (who extracted the binary algorithm from Donald Knuth,
 *  "The Art of Computer Programming", pp.441-462, vol 2, 1969).
 * 
 * Forms the power as    p = prod b_i x^(2^i), where b_i are the bits in n
 *
 * Of course, you are in trouble if the result doesn't fit in a long
 */

long pow( long x, long n )
{
    if( x==0 || n<0 ) {
        return 0;
    }

    long p = (n%2) ? x : 1;        // i=0 term
    while( n >>= 1 ) {             // Do terms i=2,...
        x *= x;                      // x^(2^i)
        if (n%2) {p *= x;}
    }

    return p;
}


double pow( double x, long n )
{
    if (x==0) return 0;                         

    if (n<0) return pow(1.0/x,-n);  

    double p = (n%2) ? x : 1;      // i=0 term
    while( n >>= 1 ) {             // Do terms i=2,...
        x *= x;                       // x^(2^i)
        if (n%2) {p *= x;}
    }

    return p;
}

/*
 * Variable number of arguments min and max
 *
 * All these functions have an extra dummy argument because they
 * are done using varargs in Cobalt Blue's run time library.  varargs
 * scares me, so I just created the required explicit functions.
 */

long   vmax(long x1, long x2, long x3, long)
{
    long q = max(x1,x2);
    return max(q,x3);
}                          

long   vmax(long x1, long x2, long x3, long x4, long)
{     
    long q = vmax(x1,x2,x3,0);
    return max(q,x4);
}

long   vmax(long x1, long x2, long x3, long x4, long x5, long)
{
    long q = vmax(x1,x2,x3,x4,0);
    return max(q,x5);
}

double   vmax(double x1, double x2, double x3, double)
{
    double q = max(x1,x2);
    return max(q,x3);
}                          

double   vmax(double x1, double x2, double x3, double x4, double)
{     
    double q = vmax(x1,x2,x3,0);
    return max(q,x4);
}

double   vmax(double x1, double x2, double x3, double x4, double x5, double)
{
    double q = vmax(x1,x2,x3,x4,0);
    return max(q,x5);
}

long   vmin(long x1, long x2, long x3, long)
{
    long q = min(x1,x2);
    return min(q,x3);
}                          

long   vmin(long x1, long x2, long x3, long x4, long)
{     
    long q = vmin(x1,x2,x3,0);
    return min(q,x4);
}

long   vmin(long x1, long x2, long x3, long x4, long x5, long)
{
    long q = vmin(x1,x2,x3,x4,0);
    return min(q,x5);
}

double   vmin(double x1, double x2, double x3, double)
{
    double q = min(x1,x2);
    return min(q,x3);
}                          

double   vmin(double x1, double x2, double x3, double x4, double)
{     
    double q = vmin(x1,x2,x3,0);
    return min(q,x4);
}

double   vmin(double x1, double x2, double x3, double x4, double x5, double)
{
    double q = vmin(x1,x2,x3,x4,0);
    return min(q,x5);
}

/*
 * Some minor string support
 *
 * f_concat concats the strings s1 and s2 into s.  The long is a dummy;
 * it is needed by Cobalt Blue because the Cobalt Blue version of f_concat
 * uses varargs.
 */

void ini_chrtmp(CHRTMP *ct, int n)  // initialize the chrtmp array 
{
    for(int i=0; i<n; i++) {
        ct[i].siz = 0;
        ct[i].s = 0; 
    }
}

void rel_chrtmp(CHRTMP *ct, int n)  // release space used by chrtmp array
{
    for(int i=0; i<n; i++) {
        delete [] ct[i].s;
        ct[i].siz = 0;
    }
}

char *f_concat(CHRTMP *target, char *s1, char *s2, long)
{
    int l1 = strlen(s1);
    int l2 = strlen(s2);
    int l = l1+l2;
    if (target->siz<l) {
        delete target->s;
        target->siz = l;
        target->s = new char [l+1];
    }
    strcpy(target->s,s1);
    strcpy(target->s+l1,s2);
    return target->s;
}
