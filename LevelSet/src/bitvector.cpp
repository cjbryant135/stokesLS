#include <stdlib.h>
#include "bitvector.h"

namespace levelset {

    BitVector::BitVector(const unsigned int m) : data(NULL)
    {
        Resize(m);
        for (int i=0; i<bytesize; ++i) mask[i] = 0x01 << i;
    }

    void BitVector::SetBit(const unsigned int i, const char b)
    {
        data[i/bytesize] = b ? data[i/bytesize] | mask[i%bytesize]
            : data[i/bytesize] & ~mask[i%bytesize];
    }

    char BitVector::ReadBit(const unsigned int i) const
    {
        return data[i/bytesize] & mask[i%bytesize];
    }

    void BitVector::Resize(const unsigned int m)
    {
        if (data) delete[] data;   
        datasize = (m+bytesize-1)/bytesize;
        data = new unsigned char[datasize];
        maxi = m;
        Clear();
    }

    std::ostream& operator<<(std::ostream& s, const BitVector& v)
    {
        s << '[';
        for (int i=0; i<v.maxi; ++i)
            s << (v.ReadBit(i) ? '1' : '0');
        s << "]\n";
        return s;
    }

}
