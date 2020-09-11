//
//  pm2io.cpp
//  LevelSet
//
//  Created by David Chopp on 7/24/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "defs.h"
#include "polarmesh2d.h"
#include "pm2boundary.h"

namespace levelset {
    
    void PolarMesh2D::WriteBinaryXY(const std::string fname, const int kstart,
                                  const int kend, const char* comment)
    {
        WriteBinaryXY(fname.c_str(), kstart, kend, comment);
    }
    
    void PolarMesh2D::WriteBinaryXY(const char* fname, const int kstart,
                                  const int kend, const char* comment)
    {
        int kval[2];
        kval[0] = kstart;
        kval[1] = kend;
        if (kval[0] == -1) {
            kval[0] = 0;
            kval[1] = maxk-1;
        }
        if (kval[1] == -1) kval[1] = kval[0];
        
        std::ofstream out(fname, std::ios::out | std::ios::binary);
        out.write(reinterpret_cast<char*>(&tmaxi),sizeof(tmaxi));
        out.write(reinterpret_cast<char*>(&tmaxj),sizeof(tmaxj));
        int klen = kval[1]-kval[0]+1;
        out.write(reinterpret_cast<char*>(&klen),sizeof(klen));
        out.write(reinterpret_cast<char*>(&(thedata[kval[0]*tmaxi*tmaxj])),
                  klen*tmaxi*tmaxj*sizeof(double));
        double* coords = new double[2*tmaxi*tmaxj];
        int xlo = bc->Width(BDRY2D_XLO);
        int ylo = bc->Width(BDRY2D_YLO);
        for (int i=0; i<tmaxi; ++i)
            for (int j=0; j<tmaxj; ++j)
                coords[i*tmaxj+j] = X(i-xlo, j-ylo);
        for (int i=0; i<tmaxi; ++i)
            for (int j=0; j<tmaxj; ++j)
                coords[tmaxi*tmaxj+i*tmaxj+j] = Y(i-xlo, j-ylo);
        out.write(reinterpret_cast<char*>(coords),(2*tmaxi*tmaxj)*sizeof(double));
        out.close();
        delete[] coords;
    }
    
    void PolarMesh2D::ReadBinaryXY(const std::string fname, const int kstart)
    {
        ReadBinaryXY(fname.c_str(), kstart);
    }
    
    void PolarMesh2D::ReadBinaryXY(const char* fname, const int kstart)
    {
        std::ifstream in(fname, std::ios::in | std::ios::binary);
        int n;
        in.read(reinterpret_cast<char*>(&n),sizeof(n));
        if (n == tmaxi) {
            in.read(reinterpret_cast<char*>(&n),sizeof(n));
            if (n == tmaxj) {
                in.read(reinterpret_cast<char*>(&n),sizeof(n));
                in.read(reinterpret_cast<char*>(&(thedata[kstart*tmaxi*tmaxj])),
                        n*tmaxi*tmaxj*sizeof(double));
                in.close();
            } else {
                std::cerr << "Error reading binary file: tmaxj doesn't match\n";
            }
        } else {
            std::cerr << "Error reading binary file: tmaxi doesn't match\n";
            exit(1);
        }
        
    }
    
    void PolarMesh2D::WriteMatlabCoordinatesXY(const char* fname)
    {
        std::ofstream out(fname, std::ios::out | std::ios::binary);
        out.write(reinterpret_cast<char*>(&tmaxi),sizeof(tmaxi));
        out.write(reinterpret_cast<char*>(&tmaxj),sizeof(tmaxj));
        
        double* coords = new double[2*tmaxi*tmaxj];
        int xlo = bc->Width(BDRY2D_XLO);
        int ylo = bc->Width(BDRY2D_YLO);
        for (int i=0; i<tmaxi; ++i)
            for (int j=0; j<tmaxj; ++j)
                coords[i*tmaxj+j] = X(i-xlo, j-ylo);
        for (int i=0; i<tmaxi; ++i)
            for (int j=0; j<tmaxj; ++j)
                coords[tmaxi*tmaxj+i*tmaxj+j] = Y(i-xlo, j-ylo);
        out.write(reinterpret_cast<char*>(coords),(2*tmaxi*tmaxj)*sizeof(double));
        out.close();
        delete[] coords;
    }
    
    void PolarMesh2D::WriteBinaryRT(const std::string fname, const int kstart,
                                    const int kend, const char* comment)
    {
        WriteBinaryRT(fname.c_str(), kstart, kend, comment);
    }
    
    void PolarMesh2D::WriteBinaryRT(const char* fname, const int kstart,
                                    const int kend, const char* comment)
    {
        int kval[2];
        kval[0] = kstart;
        kval[1] = kend;
        if (kval[0] == -1) {
            kval[0] = 0;
            kval[1] = maxk-1;
        }
        if (kval[1] == -1) kval[1] = kval[0];
        
        std::ofstream out(fname, std::ios::out | std::ios::binary);
        out.write(reinterpret_cast<char*>(&tmaxi),sizeof(tmaxi));
        out.write(reinterpret_cast<char*>(&tmaxj),sizeof(tmaxj));
        int klen = kval[1]-kval[0]+1;
        out.write(reinterpret_cast<char*>(&klen),sizeof(klen));
        out.write(reinterpret_cast<char*>(&(thedata[kval[0]*tmaxi*tmaxj])),
                  klen*tmaxi*tmaxj*sizeof(double));
        out.write(reinterpret_cast<char*>(r),tmaxi*sizeof(double));
        double* coords = new double[tmaxj];
        int ylo = bc->Width(BDRY2D_YLO);
        for (int j=0; j<tmaxj; ++j)
            coords[j] = Theta(j-ylo);
        out.write(reinterpret_cast<char*>(coords),tmaxj*sizeof(double));
        out.close();
        delete[] coords;
    }
    
    void PolarMesh2D::ReadBinaryRT(const std::string fname, const int kstart)
    {
        ReadBinaryRT(fname.c_str(), kstart);
    }
    
    void PolarMesh2D::ReadBinaryRT(const char* fname, const int kstart)
    {
        std::ifstream in(fname, std::ios::in | std::ios::binary);
        int n;
        in.read(reinterpret_cast<char*>(&n),sizeof(n));
        if (n == tmaxi) {
            in.read(reinterpret_cast<char*>(&n),sizeof(n));
            if (n == tmaxj) {
                in.read(reinterpret_cast<char*>(&n),sizeof(n));
                in.read(reinterpret_cast<char*>(&(thedata[kstart*tmaxi*tmaxj])),
                        n*tmaxi*tmaxj*sizeof(double));
                in.close();
            } else {
                std::cerr << "Error reading binary file: tmaxj doesn't match\n";
            }
        } else {
            std::cerr << "Error reading binary file: tmaxi doesn't match\n";
            exit(1);
        }
        
    }
    
    void PolarMesh2D::WriteMatlabCoordinatesRT(const char* fname)
    {
        std::ofstream out(fname, std::ios::out | std::ios::binary);
        out.write(reinterpret_cast<char*>(&tmaxi),sizeof(tmaxi));
        out.write(reinterpret_cast<char*>(&tmaxj),sizeof(tmaxj));
        
        out.write(reinterpret_cast<char*>(r),tmaxi*sizeof(double));
        double* coords = new double[tmaxj];
        int ylo = bc->Width(BDRY2D_YLO);
        for (int j=0; j<tmaxj; ++j)
            coords[j] = Theta(j-ylo);
        out.write(reinterpret_cast<char*>(coords),tmaxj*sizeof(double));
        out.close();
        delete[] coords;
    }

}
