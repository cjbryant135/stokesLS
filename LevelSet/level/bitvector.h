#ifndef __BIT_VECTOR_H__
#define __BIT_VECTOR_H__

#include <iostream>

namespace levelset {

#ifndef BYTESIZE
#define BYTESIZE
    const unsigned short bytesize = sizeof(unsigned char)*8;
#endif

    class BitVector 
    {
    private:
   
        unsigned char *data;
        unsigned int  datasize;
        unsigned char mask[bytesize];
        unsigned int  maxi;

    public:

        BitVector(const unsigned int m);
        ~BitVector(void) {if (data) delete[] data;}

        void Clear(void) {for (int i=0; i<datasize; ++i) data[i] = 0x00;} 
        void SetBit(const unsigned int i, const char b = 0x01);
        char ReadBit(const unsigned int i) const;

        void Resize(const unsigned int m);

        friend std::ostream& operator<<(std::ostream& s, const BitVector& v);
//   friend istream& operator>>(istream& s, BitVector& v);
   
    };

}

#endif
