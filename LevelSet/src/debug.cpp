#include "debug.h"
#include <iostream>
#include <iomanip>

namespace levelset {

    void ShowStatus(std::ostream& s, const int* array, const int mi, const int mj,
                    const int k)
    {
        for (int j=mj-1; j>=0; --j) {
            for (int i=0; i<mi; ++i)
                if (array[k*mi*mj+i*mj+j])
                    s << (char)('A' + array[k*mi*mj+i*mj+j]);
                else
                    s << ' ';
            s << '\n';
        }
    }

    void PlusMinus(std::ostream& s, const double* array, const int mi, const int mj)
    {
        for (int j=mj-1; j>=0; --j) {
            for (int i=0; i<mi; ++i)
                s << (array[mj*i+j] > 0 ? '+' : (array[mj*i+j] < 0 ? '-' : '0'));
            s << '\n';
        }
    }

#ifdef USE_MATHERR

    int matherr(struct exception *x)
    {
        cerr << "Math exception detected:\n";
        switch (x->type) {
        case DOMAIN:
            cerr << "Domain error in ";
            break;
        case SING:
            cerr << "Singularity in ";
            break;
        case OVERFLOW:
            cerr << "Overflow in ";
            break;
        case UNDERFLOW:
            cerr << "Underflow in ";
            break;
        }
        cerr << x->name << "\targ1 = " << x->arg1
             << " arg2 = " << x->arg2
             << " retval = " << x->retval;
        return 1;
    }

#endif

}
