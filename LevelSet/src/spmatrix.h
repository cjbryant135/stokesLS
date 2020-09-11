#ifndef __SPMATRIX_H__
#define __SPMATRIX_H__

#include "spvector.h"
#include "eltype.h"
#include <iostream>

namespace blockmatrix {
        
    class spMatrix
    {
    public:

        size_t RowDim;
        size_t ClmDim;
        size_t *Len;
        ElType **El;
        myBool ElSorted;

        spMatrix(const size_t r, const size_t c);
        virtual ~spMatrix(void);
   
        virtual void Clear(void);

        virtual void reshape(const size_t r, const size_t c);
        size_t rows(void) const {return RowDim;} 
        size_t cols(void) const {return ClmDim;} 

        double Value(const size_t r, const size_t c) const;
        inline double operator()(const size_t r, const size_t c) const {return Value(r,c);}

        double& operator()(const size_t r, const size_t c);

        void Zero(const size_t r, const size_t c);

        void ExtendDim(const size_t r, const size_t c);
        void ExtendRow(const size_t r, const size_t e);
        void Allocate(const size_t r, const size_t e);

        double* ColSums(double* cols = NULL);

        spMatrix& operator=(const spMatrix& A);

        void print(std::ostream& s = std::cout) const;
    };

    spVector operator*(const spMatrix& A, const spVector& X);

}
#endif


