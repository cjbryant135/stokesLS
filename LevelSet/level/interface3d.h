#ifndef __INTERFACE3D_H__
#define __INTERFACE3D_H__

#include "datavec.h"
#ifndef NO_GRAPHICS
#include "plotwindow3d.h"
#endif
#include "segment.h"
#include "initialfunc.h"
#include "uniformmesh3d.h"

namespace levelset {
        
    class Node {
    public:

        double x, y, z;
        
    Node(const double a, const double b, const double c) 
        : x(a), y(b), z(c) {}
    Node(void) : x(0.), y(0.), z(0.) {}
        
        void Set(const double a, const double b, const double c)
        {x=a; y=b; z=c;}
    };

    class Element {
    public:

        int node[3];
        Element(const int i, const int j, const int k) {Set(i,j,k);}
        Element(void) {Set(-1,-1,-1);}

        void Set(const int i, const int j, const int k)
        {node[0] = i; node[1] = j, node[2] = k;}
        
        void Set(const int* n)
        {for (int i=0; i<3; ++i) node[i] = n[i];}
    };

    class UniformMesh3D;

    class Interface3D {
    public:

        DataVec<Node> node;
        DataVec<Element> elem;
        int ncount, ecount;

    Interface3D(void) : node(), elem(), ncount(0), ecount(0) {;} 

        void Trim(const UniformMesh3D& m, const int ksurf, const int kfilt, 
                  const int sval = 0, const double fval = 0, const int sgn = 1, 
                  const double thresh = 1.0e-4, const char adjust = 1);
        
        void Dump(const char* filename) const;

#ifndef NO_GRAPHICS
        void Plot(PlotWindow3D& port, const char fill = 0) const ;
#endif
    };

}

#endif // __INTERFACE3D_H__

