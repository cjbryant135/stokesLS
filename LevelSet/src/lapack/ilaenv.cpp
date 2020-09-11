/*
 * C++ implementation of lapack routine ilaenv
 *
 * $Id: ilaenv.cpp,v 1.6 1993/07/08 22:35:09 alv Exp $
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
 * Translated output further modified to remove Fortran string
 * processing.
 *
 * Translated by FOR_C++, v1.1 (P), on 02/18/93 at 07:27:37
 * FOR_C++ Options SET: alloc do=rt no=p pf=xlapack s=dv str=l - prototypes
 *
 * $Log: ilaenv.cpp,v $
 * Revision 1.6  1993/07/08 22:35:09  alv
 * recognizes RW_USE_UNBLOCKED
 *
 * Revision 1.5  1993/04/06  20:42:56  alv
 * added const to parameters; added include lapkdefs
 *
 * Revision 1.4  1993/03/19  18:47:27  alv
 * commented out unused vars to shut up SUN warnings
 *
 * Revision 1.3  1993/03/19  17:18:24  alv
 * added RWLAPKDECL linkage specifier
 *
 * Revision 1.2  1993/03/05  23:18:01  alv
 * changed ref parms to const ref
 *
 * Revision 1.1  1993/03/03  16:09:38  alv
 * Initial revision
 *
 */

#include <ctype.h>
#include <string.h>
#if 0
#include "rw/lapkdefs.h"
#include "rw/bla.h"
#include "rw/lapack.h"
#include "rw/fortran.h" /* Fortran run time library */
#else
#include "level/lapack.h"
#endif

RWLAPKDECL long /*FUNCTION*/ ilaenv(const long &ispec, char *name, char * /*opts*/, 
                                    const long &n1, const long &n2, const long &/*n3*/, const long &n4)
{
    char c2[3], c3[4], c4[3], subnam[7];
    int cname, sname;
    char c1;
    long /*i, i_, ic,*/ ilaenv_v, /*iz,*/ nb, nbmin, nx;

  
    //  -- LAPACK auxiliary routine (preliminary version) --
    //     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    //     Courant Institute, Argonne National Lab, and Rice University
    //     February 20, 1992
  
    //     .. Scalar Arguments ..
    //     ..
  
    //  Purpose
    //  =======
  
    //  ILAENV is called from the LAPACK routines to choose problem-dependent
    //  parameters for the local environment.  See ISPEC for a description of
    //  the parameters.
  
    //  This version provides a set of parameters which should give good,
    //  but not optimal, performance on many of the currently available
    //  computers.  Users are encouraged to modify this subroutine to set
    //  the tuning parameters for their particular machine using the option
    //  and problem size information in the arguments.
  
    //  This routine will not function correctly if it is converted to all
    //  lower case.  Converting it to all upper case is allowed.
  
    //  Arguments
    //  =========
  
    //  ISPEC   (input) INTEGER
    //          Specifies the parameter to be returned as the value of
    //          ILAENV.
    //          = 1: the optimal blocksize; if this value is 1, an unblocked
    //               algorithm will give the best performance.
    //          = 2: the minimum block size for which the block routine
    //               should be used; if the usable block size is less than
    //               this value, an unblocked routine should be used.
    //          = 3: the crossover point (in a block routine, for N less
    //               than this value, an unblocked routine should be used)
    //          = 4: the number of shifts, used in the nonsymmetric
    //               eigenvalue routines
    //          = 5: the minimum column dimension for blocking to be used;
    //               rectangular blocks must have dimension at least k by m,
    //               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
    //          = 6: the crossover point for the SVD (when reducing an m by n
    //               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
    //               this value, a QR factorization is used first to reduce
    //               the matrix to a triangular form.)
    //          = 7: the number of processors
    //          = 8: the crossover point for the multishift QR and QZ methods
    //               for nonsymmetric eigenvalue problems.
  
    //  NAME    (input) CHARACTER*(*)
    //          The name of the calling subroutine, in either upper case or
    //          lower case.
  
    //  OPTS    (input) CHARACTER*(*)
    //          The character options to the subroutine NAME, concatenated
    //          into a single character string.  For example, UPLO = 'U',
    //          TRANS = 'T', and DIAG = 'N' for a triangular routine would
    //          be specified as OPTS = 'UTN'.
  
    //  N1      (input) INTEGER
    //  N2      (input) INTEGER
    //  N3      (input) INTEGER
    //  N4      (input) INTEGER
    //          Problem dimensions for the subroutine NAME; these may not all
    //          be required.
  
    // (ILAENV) (output) INTEGER
    //          >= 0: the value of the parameter specified by ISPEC
    //          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
  
    //  Further Details
    //  ===============
  
    //  The following conventions have been used when calling ILAENV from the
    //  LAPACK routines:
    //  1)  OPTS is a concatenation of all of the character options to
    //      subroutine NAME, in the same order that they appear in the
    //      argument list for NAME, even if they are not used in determining
    //      the value of the parameter specified by ISPEC.
    //  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
    //      that they appear in the argument list for NAME.  N1 is used
    //      first, N2 second, and so on, and unused problem dimensions are
    //      passed a value of -1.
    //  3)  The parameter value returned by ILAENV is checked for validity in
    //      the calling subroutine.  For example, ILAENV is used to retrieve
    //      the optimal blocksize for STRTRI as follows:
  
    //      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
    //      IF( NB.LE.1 ) NB = MAX( 1, N )
  
    //  =====================================================================
  
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..

#ifdef RW_USE_UNBLOCKED
    if (ispec==1) return 1;   // optimal block size always one
#endif
  
    switch( ispec ) { 
    case 1: goto L_100;
    case 2: goto L_100;
    case 3: goto L_100;
    case 4: goto L_400;
    case 5: goto L_500;
    case 6: goto L_600;
    case 7: goto L_700;
    case 8: goto L_800;
    }
  
    //     Invalid value for ISPEC
  
    ilaenv_v = -1;
    return( ilaenv_v );
  
L_100:
    ;
  
    //     Convert NAME to upper case if the first character is lower case.
    //     Store upper case name in subnam.

    ilaenv_v = 1;
    subnam[0] = toupper(name[0]);
    subnam[1] = toupper(name[1]);
    subnam[2] = toupper(name[2]);
    subnam[3] = toupper(name[3]);
    subnam[4] = toupper(name[4]);
    subnam[5] = toupper(name[5]);
    subnam[6] = 0;
  
    c1 = subnam[0];
    sname = c1 == 'S' || c1 == 'D';
    cname = c1 == 'C' || c1 == 'Z';
    if( !(cname || sname) ) 
        return( ilaenv_v );

    c2[0] = subnam[1];    // Fortran: c2 = subnam(2:3)
    c2[1] = subnam[2];
    c2[2] = 0;
    c3[0] = subnam[3];    // Fortran: c3 = subnam(4:6)
    c3[1] = subnam[4];
    c3[2] = subnam[5];
    c3[3] = 0;
    c4[0] = c3[1];        // Fortran: c4 = c3(2:3)
    c4[1] = c3[2];
    c4[2] = 0;
  
    switch( ispec ) { 
    case 1: goto L_110;
    case 2: goto L_200;
    case 3: goto L_300;
    }
  
L_110:
    ;
  
    //     ISPEC = 1:  block size
  
    //     In these examples, separate code is provided for setting NB for
    //     real and DComplex.  We assume that NB will take the same value in
    //     single or double precision.
  
    nb = 1;
  
    if( strcmp(c2,"GE") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            if( sname ) { 
                nb = 64;
            }
            else { 
                nb = 64;
            }
        }
        else if( ((strcmp(c3,"QRF") == 0 || strcmp(c3,"RQF") == 0) || 
                  strcmp(c3,"LQF") == 0) || strcmp(c3,"QLF") == 0 ) { 
            if( sname ) { 
                nb = 32;
            }
            else { 
                nb = 32;
            }
        }
        else if( strcmp(c3,"HRD") == 0 ) { 
            if( sname ) { 
                nb = 32;
            }
            else { 
                nb = 32;
            }
        }
        else if( strcmp(c3,"BRD") == 0 ) { 
            if( sname ) { 
                nb = 32;
            }
            else { 
                nb = 32;
            }
        }
        else if( strcmp(c3,"TRI") == 0 ) { 
            if( sname ) { 
                nb = 64;
            }
            else { 
                nb = 64;
            }
        }
    }
    else if( strcmp(c2,"PO") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            if( sname ) { 
                nb = 64;
            }
            else { 
                nb = 64;
            }
        }
    }
    else if( strcmp(c2,"SY") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            if( sname ) { 
                nb = 64;
            }
            else { 
                nb = 64;
            }
        }
        else if( sname && strcmp(c3,"TRD") == 0 ) { 
            nb = 1;
        }
        else if( sname && strcmp(c3,"GST") == 0 ) { 
            nb = 64;
        }
    }
    else if( cname && strcmp(c2,"HE") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            nb = 64;
        }
        else if( strcmp(c3,"TRD") == 0 ) { 
            nb = 1;
        }
        else if( strcmp(c3,"GST") == 0 ) { 
            nb = 64;
        }
    }
    else if( sname && strcmp(c2,"OR") == 0 ) { 
        if( c3[0] == 'G' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nb = 32;
            }
        }
        else if( c3[0] == 'M' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nb = 32;
            }
        }
    }
    else if( cname && strcmp(c2,"UN") == 0 ) { 
        if( c3[0] == 'G' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nb = 32;
            }
        }
        else if( c3[0] == 'M' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nb = 32;
            }
        }
    }
    else if( strcmp(c2,"GB") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            if( sname ) { 
                if( n4 <= 64 ) { 
                    nb = 1;
                }
                else { 
                    nb = 32;
                }
            }
            else { 
                if( n4 <= 64 ) { 
                    nb = 1;
                }
                else { 
                    nb = 32;
                }
            }
        }
    }
    else if( strcmp(c2,"PB") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            if( sname ) { 
                if( n2 <= 64 ) { 
                    nb = 1;
                }
                else { 
                    nb = 32;
                }
            }
            else { 
                if( n2 <= 64 ) { 
                    nb = 1;
                }
                else { 
                    nb = 32;
                }
            }
        }
    }
    else if( strcmp(c2,"TR") == 0 ) { 
        if( strcmp(c3,"TRI") == 0 ) { 
            if( sname ) { 
                nb = 64;
            }
            else { 
                nb = 64;
            }
        }
    }
    else if( strcmp(c2,"LA") == 0 ) { 
        if( strcmp(c3,"UUM") == 0 ) { 
            if( sname ) { 
                nb = 64;
            }
            else { 
                nb = 64;
            }
        }
    }
    else if( sname && strcmp(c2,"ST") == 0 ) { 
        if( strcmp(c3,"EBZ") == 0 ) { 
            nb = 1;
        }
    }
    ilaenv_v = nb;
    return( ilaenv_v );
  
L_200:
    ;
  
    //     ISPEC = 2:  minimum block size
  
    nbmin = 2;
    if( strcmp(c2,"GE") == 0 ) { 
        if( ((strcmp(c3,"QRF") == 0 || strcmp(c3,"RQF") == 0) || strcmp(c3
                                                                        ,"LQF") == 0) || strcmp(c3,"QLF") == 0 ) { 
            if( sname ) { 
                nbmin = 2;
            }
            else { 
                nbmin = 2;
            }
        }
        else if( strcmp(c3,"HRD") == 0 ) { 
            if( sname ) { 
                nbmin = 2;
            }
            else { 
                nbmin = 2;
            }
        }
        else if( strcmp(c3,"BRD") == 0 ) { 
            if( sname ) { 
                nbmin = 2;
            }
            else { 
                nbmin = 2;
            }
        }
        else if( strcmp(c3,"TRI") == 0 ) { 
            if( sname ) { 
                nbmin = 2;
            }
            else { 
                nbmin = 2;
            }
        }
    }
    else if( strcmp(c2,"SY") == 0 ) { 
        if( strcmp(c3,"TRF") == 0 ) { 
            if( sname ) { 
                nbmin = 2;
            }
            else { 
                nbmin = 2;
            }
        }
        else if( sname && strcmp(c3,"TRD") == 0 ) { 
            nbmin = 2;
        }
    }
    else if( cname && strcmp(c2,"HE") == 0 ) { 
        if( strcmp(c3,"TRD") == 0 ) { 
            nbmin = 2;
        }
    }
    else if( sname && strcmp(c2,"OR") == 0 ) { 
        if( c3[0] == 'G' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nbmin = 2;
            }
        }
        else if( c3[0] == 'M' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nbmin = 2;
            }
        }
    }
    else if( cname && strcmp(c2,"UN") == 0 ) { 
        if( c3[0] == 'G' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nbmin = 2;
            }
        }
        else if( c3[0] == 'M' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nbmin = 2;
            }
        }
    }
    ilaenv_v = nbmin;
    return( ilaenv_v );
  
L_300:
    ;
  
    //     ISPEC = 3:  crossover point
  
    nx = 0;
    if( strcmp(c2,"GE") == 0 ) { 
        if( ((strcmp(c3,"QRF") == 0 || strcmp(c3,"RQF") == 0) || strcmp(c3
                                                                        ,"LQF") == 0) || strcmp(c3,"QLF") == 0 ) { 
            if( sname ) { 
                nx = 128;
            }
            else { 
                nx = 128;
            }
        }
        else if( strcmp(c3,"HRD") == 0 ) { 
            if( sname ) { 
                nx = 128;
            }
            else { 
                nx = 128;
            }
        }
        else if( strcmp(c3,"BRD") == 0 ) { 
            if( sname ) { 
                nx = 128;
            }
            else { 
                nx = 128;
            }
        }
    }
    else if( strcmp(c2,"SY") == 0 ) { 
        if( sname && strcmp(c3,"TRD") == 0 ) { 
            nx = 1;
        }
    }
    else if( cname && strcmp(c2,"HE") == 0 ) { 
        if( strcmp(c3,"TRD") == 0 ) { 
            nx = 1;
        }
    }
    else if( sname && strcmp(c2,"OR") == 0 ) { 
        if( c3[0] == 'G' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nx = 128;
            }
        }
    }
    else if( cname && strcmp(c2,"UN") == 0 ) { 
        if( c3[0] == 'G' ) { 
            if( (((((strcmp(c4,"QR") == 0 || strcmp(c4,"RQ") == 0) || 
                    strcmp(c4,"LQ") == 0) || strcmp(c4,"QL") == 0) || strcmp(c4
                                                                             ,"HR") == 0) || strcmp(c4,"TR") == 0) || strcmp(c4,"BR") == 0
                ) { 
                nx = 128;
            }
        }
    }
    ilaenv_v = nx;
    return( ilaenv_v );
  
L_400:
    ;
  
    //     ISPEC = 4:  number of shifts (used by xHSEQR)
  
    ilaenv_v = 6;
    return( ilaenv_v );
  
L_500:
    ;
  
    //     ISPEC = 5:  minimum column dimension (not used)
  
    ilaenv_v = 2;
    return( ilaenv_v );
  
L_600:
    ;
  
    //     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
  
    ilaenv_v = (long)( (float)( min( n1, n2 ) )*1.6e0 );
    return( ilaenv_v );
  
L_700:
    ;
  
    //     ISPEC = 7:  number of processors (not used)
  
    ilaenv_v = 1;
    return( ilaenv_v );
  
L_800:
    ;
  
    //     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
  
    ilaenv_v = 50;
    return( ilaenv_v );
  
    //     End of ILAENV
  
} // end of function 

