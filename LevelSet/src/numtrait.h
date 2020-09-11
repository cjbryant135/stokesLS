#ifndef __NUMTRAIT_H__
#define __NUMTRAIT_H__

#include <math.h>
#include <limits.h>
#include <float.h>

namespace levelset {
        
    template <class T>
        class NumericTrait {
    public:
   
        const T max_value;
        const T min_value;

    NumericTrait(void) : max_value(INT_MAX), min_value(-INT_MAX) {;}
    };

    template<>
        class NumericTrait<double> {
    public:
   
        const double max_value;
        const double min_value;

    NumericTrait(void) : max_value(DBL_MAX), min_value(-DBL_MAX) {;}
    };

#if 0
    class NumericTrait<float> {
    public:
   
        const float max_value = FLT_MAX;
        const float min_value = -FLT_MAX;
        const char tname[20] = "float";
    };

    class NumericTrait<int> {
    public:
   
        const int max_value = INT_MAX;
        const int min_value = -INT_MAX;
        const char tname[20] = "int";
    };
#endif

}
#endif // __NUMTRAIT_H__
