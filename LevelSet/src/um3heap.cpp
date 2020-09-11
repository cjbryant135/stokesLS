#include <iostream>
#include <iomanip>
#include "um3heap.h"
#ifdef MEMWATCH
#include "memwatch.h"
#endif

namespace levelset {
        
    void UniformMesh3D_Heap::Change(const int i, const double val, const int ikhi)
    {
        data[i] = val;
        SortFrom(i,ikhi);
    }

    void UniformMesh3D_Heap::Insert(const int i, const int j, const int k,
                                    const double val, const int ikhi)
    {
        mesh->idata_(i,j,k,ikhi) = length;
        ++length;
        if (length > data.Length()) {
            data.Resize(heapT_row2(length));
        }
        data[length-1] = HeapElement3D(i,j,k,val);
        SortFrom(length-1,ikhi);
    }

    void UniformMesh3D_Heap::Extract(int& i, int& j, int& k, const int ikhi)
    {
        HeapElement3D he = data[0];
        --length;
        if (length > 0) {
            MoveElt(length,0,ikhi);
            SortFrom(0,ikhi);
        }
        i = he.Index(0);
        j = he.Index(1);
        k = he.Index(2);
    }

    void UniformMesh3D_Heap::SortFrom(unsigned int n, const int ikhi)
    {
        if (type == MinHeap) {
            int k;
            for (k=n; k>0; k=parent(k)) 
                if (data[k] < data[parent(k)]) 
                    SwapElts(k,parent(k),ikhi);
                else 
                    break;

            int c;
            while (child(k,HLeft) < length) {
                if (child(k,HRight) >= length) 
                    c = child(k,HLeft);
                else
                    c = (data[child(k,HLeft)] < data[child(k,HRight)]) ? 
                        child(k,HLeft) : child(k,HRight);
                if (data[k] > data[c]) {
                    SwapElts(k,c,ikhi);
                    k = c;
                } else 
                    break;
            }
        } else {
            int k;
            for (k=n; k>0; k=parent(k)) 
                if (data[k] > data[parent(k)]) 
                    SwapElts(k,parent(k),ikhi);
                else 
                    break;

            int c;
            while (child(k,HLeft) < length) {
                if (child(k,HRight) >= length) 
                    c = child(k,HLeft);
                else
                    c = (data[child(k,HLeft)] > data[child(k,HRight)]) ? 
                        child(k,HLeft) : child(k,HRight);
                if (data[k] < data[c]) {
                    SwapElts(k,c,ikhi);
                    k = c;
                } else 
                    break;
            }
        }

    }

#ifdef LEVEL_DEBUG
    void UniformMesh3D_Heap::Dump(std::ostream& s, const char val) const
    {
        int depth = heapT_log2(length);
        int rows = 0x01 << (depth-1);
        int width = 9;
        int prec = 3;

        s << std::setfill(' ');
        for (int row=0; row<rows; ++row) {
            for (int d=0; d<depth; ++d) 
                if (!(row%(0x01 << (depth-d-1)))) {
                    int b = (0x01 << d) - 1;
                    int r = row/(0x01 << (depth-d-1));
                    if (b+r < length) {
                        if (val) {
                            s << std::setw(width) << std::setprecision(prec) << data[b+r].Value();
                        }
                        else {
                            s << std::setw(width/3) << data[b+r].Index(0) << ',' \
                              << std::setw(width/3) << data[b+r].Index(1) << ',' \
                              << std::setw(width/3) << data[b+r].Index(2)        \
                              << std::setw(width-2*(width/3)-1) << " ";
                        }
                    }
                    else
                        s << std::setw(width) << " ";
                } else 
                    s << std::setw(width) << " ";
            s << "\n";
        }
    }

    void UniformMesh3D_Heap::Check(const int ikhi) const
    {
        for (int ii=0; ii<length; ++ii) 
            if (ii != mesh->idata_(data[ii].Index(0),data[ii].Index(1),
                                   data[ii].Index(2),ikhi)) {
                std::cout << "Error in um3heap.cpp: \n";
                std::cout << "Heap Index = " << ii << ", and idata_(" << data[ii].Index(0);
                std::cout << "," << data[ii].Index(1) << ",";
                std::cout << "," << data[ii].Index(2) << ", ikhi) = ";
                std::cout << mesh->idata_(data[ii].Index(0),data[ii].Index(1),
                                          data[ii].Index(2),ikhi);
                std::cout << "\n";
            }
    }

#endif


}











