//#include <stdlib.h>
#include "bitmatrix.h"

namespace levelset {

    BitMatrix::BitMatrix(const unsigned int m, const unsigned int n) : data(NULL)
    {
        Resize(m, n);
        for (int i=0; i<bytesize; ++i) mask[i] = 0x01 << i;
    }

    BitMatrix::BitMatrix(void) : maxi(0), maxj(0), data(NULL)
    {
        for (int i=0; i<bytesize; ++i) mask[i] = 0x01 << i;
    }

    void BitMatrix::SetBit(const unsigned int i, const unsigned int j,
                           const char b)
    {
        int bit = maxj*i+j;
        data[bit/bytesize] = b ? data[bit/bytesize] | mask[bit%bytesize]
            : data[bit/bytesize] & ~mask[bit%bytesize];
    }

    char BitMatrix::ReadBit(const unsigned int i, const unsigned int j) const
    {
        int bit = maxj*i+j;
        return data[bit/bytesize] & mask[bit%bytesize];
    }

    void BitMatrix::Resize(const unsigned int m, const unsigned int n)
    {
        if (data) delete[] data;   
        datasize = (m*n+bytesize-1)/bytesize;
        data = new unsigned char[datasize];
        maxi = m;
        maxj = n;
        Clear();
    }

    std::ostream& operator<<(std::ostream& s, const BitMatrix& m)
    {
        for (int i=0; i<m.maxi; ++i) {
            s << '[';
            for (int j=0; j<m.maxj; ++j)
                s << (m.ReadBit(i,j) ? '1' : '0');
            s << "]\n";
        }
        return s;
    }

}
