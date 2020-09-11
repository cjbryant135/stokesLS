#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <iostream>
#include <math.h>

namespace levelset {

    void PlusMinus(std::ostream& s, const double* array, const int mi, const int mj);

    void ShowStatus(std::ostream& s, const int* array, const int mi, const int mj,
                    const int k);

#ifdef USE_MATHERR
    int matherr(struct exception *x);
#endif

#ifdef CHECK_DIVIDE_BY_ZERO
#define checkdenom(A) {if ((A)==0) std::cerr << "Divide by zero at "    \
                                             << __FILE__ << ':' << __LINE__ << '\n';}
#else
#define checkdenom(A)
#endif
        
}

#endif // __DEBUG_H__
