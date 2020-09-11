#include <iostream>
#include <iomanip>
#include "um2heap.h"

namespace levelset {
        
    void UniformMesh2D_Heap::Change(const int i, const double val, const int ikhi)
    {
        data[i] = val;
        SortFrom(i,ikhi);
    }

    void UniformMesh2D_Heap::Insert(const int i, const int j, const double val,
                                    const int ikhi)
    {
        mesh->midata_(mesh,i,j,ikhi) = length;
        ++length;
        if (length > data.Length()) {
            data.Resize(heapT_row2(length));
        }
        data[length-1] = HeapElement2D(i,j,val);
        SortFrom(length-1,ikhi);
    }

    void UniformMesh2D_Heap::Extract(int& i, int& j, const int ikhi)
    {
        HeapElement2D he = data[0];
        --length;
        if (length > 0) {
            MoveElt(length,0,ikhi);
            SortFrom(0,ikhi);
        }
        i = he.Index(0);
        j = he.Index(1);
    }

    void UniformMesh2D_Heap::SortFrom(unsigned int n, const int ikhi)
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

    void UniformMesh2D_Heap::Dump(std::ostream& s, const char val) const
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
                              << std::setw(width/3) << data[b+r].Index(1)        \
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

    void UniformMesh2D_Heap::Check(const int ikhi) const
    {
        for (int ii=0; ii<length; ++ii) 
            if (ii != mesh->midata_(mesh,data[ii].Index(0),data[ii].Index(1),ikhi)) {
                std::cout << "Error in um2heap.cpp: \n";
                std::cout << "Heap Index = " << ii << ", and idata_(" << data[ii].Index(0);
                std::cout << "," << data[ii].Index(1) << ", ikhi) = ";
                std::cout << mesh->midata_(mesh,data[ii].Index(0),data[ii].Index(1),ikhi);
                std::cout << "\n";
            }
    }

}













