#ifndef __BIT_MATRIX_H__
#define __BIT_MATRIX_H__

#include <iostream>

namespace levelset {

#ifndef BYTESIZE
#define BYTESIZE
    const unsigned short bytesize = sizeof(unsigned char)*8;
#endif

    class BitMatrix 
    {
    private:
   
        unsigned char *data;
        unsigned int  datasize;
        unsigned char mask[bytesize];
        unsigned int  maxi, maxj;

    public:

        BitMatrix(const unsigned int m, const unsigned int n);
        BitMatrix(void);
        ~BitMatrix(void) {if (data) delete[] data;}

        void Clear(void) {for (int i=0; i<datasize; ++i) data[i] = 0x00;} 
        void SetBit(const unsigned int i, const unsigned int j,
                    const char b = 0x01);
        char ReadBit(const unsigned int i, const unsigned int j) const;

        void Resize(const unsigned int m, const unsigned int n);

        friend std::ostream& operator<<(std::ostream& s, const BitMatrix& m);
    };

}

#endif
