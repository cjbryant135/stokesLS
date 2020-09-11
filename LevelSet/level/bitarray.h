#ifndef __BIT_ARRAY_H__
#define __BIT_ARRAY_H__

#include <iostream>

namespace levelset {

#ifndef BYTESIZE
#define BYTESIZE
    const unsigned short bytesize = sizeof(unsigned char)*8;
#endif

    class BitArray
    {
    private:
   
        unsigned char *data;
        unsigned int  datasize;
        unsigned char mask[bytesize];
        unsigned int  maxi, maxj, maxk;

    public:

        BitArray(const unsigned int m, const unsigned int n, const unsigned int o);
        BitArray(void);
        ~BitArray(void) {if (data) delete[] data;}

        void Clear(void) {for (int i=0; i<datasize; ++i) data[i] = 0x00;} 
        void SetBit(const unsigned int i, const unsigned int j, const unsigned int k, 
                    const char b = 0x01);
        char ReadBit(const unsigned int i, const unsigned int j, 
                     const unsigned int k) const;

        void Resize(const unsigned int m, const unsigned int n, const unsigned int o);

        friend std::ostream& operator<<(std::ostream& s, const BitArray& m);
    };

}

#endif
