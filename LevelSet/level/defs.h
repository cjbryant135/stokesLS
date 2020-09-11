#ifndef __DEFS_H__
#define __DEFS_H__

#ifndef TRUE
#define TRUE 0x01
#define FALSE 0x00
#define true 0x01
#define false 0x00
#endif

#include <math.h>

// INF - numerical representation of infinity
#ifdef __SC__
#include <SANE.h>
#define INF __inf
#endif
#ifdef __SUNPRO_CC
#include <sunmath.h>
#define INF infinity()
//#define NULL 0L
#endif

// NO_TEMPLATE_TYPEDEFS - typedefs inside template class declarations.
//                        Define this constant if templates are not allowed
//                        to contain typedefs
#ifdef __GNUC__
#define NO_TEMPLATE_TYPEDEFS 1
#endif

// NO_INCLUDE_CC_FILES - include .cc files in header files of template
//                       declarations.  Define this constant if .cc template
//                       definition files should be separate from the .h
//                       template declaration file.
//#define NO_INCLUDE_CC_FILES

// M_PI Definition
#ifndef M_PI
#define M_PI 3.1415926
#endif

//**********************************************************************
//
// End of Customizable section
//
//**********************************************************************

#ifndef NULL
#define NULL ((void*)0)
#endif

template <class T>
inline T sqr(const T x) { return x*x; }

template <class T>
inline int sign(const T x) {return x > 0 ? 1 : (x < 0 ? -1 : 0);}

#ifdef __MWERKS__
inline double abs(const double x) {return fabs(x);}
#else
#if 0
template <class T>
inline T abs(const T x) {return x >= 0 ? x : -x;}
#endif
#endif

#endif
