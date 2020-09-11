#include <stdlib.h>
#include "bitarray.h"

namespace levelset {

    BitArray::BitArray(const unsigned int m, const unsigned int n, 
                       const unsigned int o) : data(NULL)
    {
        Resize(m, n, o);
        for (int i=0; i<bytesize; ++i) mask[i] = 0x01 << i;
    }

    BitArray::BitArray(void) : maxi(0), maxj(0), maxk(0), data(NULL)
    {
        for (int i=0; i<bytesize; ++i) mask[i] = 0x01 << i;
    }

    void BitArray::SetBit(const unsigned int i, const unsigned int j,
                          const unsigned int k, const char b)
    {
        unsigned int bit = maxi*maxj*k+maxj*i+j;
        data[bit/bytesize] = b ? data[bit/bytesize] | mask[bit%bytesize]
            : data[bit/bytesize] & ~mask[bit%bytesize];
    }

    char BitArray::ReadBit(const unsigned int i, const unsigned int j,
                           const unsigned int k) const
    {
        unsigned int bit = maxi*maxj*k+maxj*i+j;
        return data[bit/bytesize] & mask[bit%bytesize];
    }

    void BitArray::Resize(const unsigned int m, const unsigned int n, 
                          const unsigned int o)
    {
        if (data) delete[] data;   
        datasize = (m*n*o+bytesize-1)/bytesize;
        data = new unsigned char[datasize];
        maxi = m;
        maxj = n;
        maxk = o;
        Clear();
    }

    std::ostream& operator<<(std::ostream& s, const BitArray& m)
    {
        for (int i=0; i<m.maxi; ++i) {
            s << '[';
            for (int j=0; j<m.maxj; ++j)  {
                s << '[';
                for (int k=0; k<m.maxk; ++k) 
                    s << (m.ReadBit(i,j,k) ? '1' : '0');
                s << "]\n";
            }
            s << "]\n";
        }
        return s;
    }

}
