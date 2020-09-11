#include "interface3d.h"
#include "numtrait.h"
#include "plotwindow3d.h"
#include "uniformmesh3d.h"
#include <sstream>
#include <fstream>

namespace levelset {
        
#ifndef NO_GRAPHICS
    void Interface3D::Plot(PlotWindow3D& port, const char fill) const
    {
        if (fill) {
            for (int i=0; i<ecount; ++i) 
                port.Triangle(node[elem[i].node[0]].x, node[elem[i].node[0]].y, node[elem[i].node[0]].z,
                              node[elem[i].node[1]].x, node[elem[i].node[1]].y, node[elem[i].node[1]].z,
                              node[elem[i].node[2]].x, node[elem[i].node[2]].y, node[elem[i].node[2]].z);
        } else {
            for (int i=0; i<ecount; ++i) {
                port.Line(node[elem[i].node[0]].x, node[elem[i].node[0]].y, node[elem[i].node[0]].z,
                          node[elem[i].node[1]].x, node[elem[i].node[1]].y, node[elem[i].node[1]].z);
                port.Line(node[elem[i].node[2]].x, node[elem[i].node[2]].y, node[elem[i].node[2]].z,
                          node[elem[i].node[1]].x, node[elem[i].node[1]].y, node[elem[i].node[1]].z);
                port.Line(node[elem[i].node[0]].x, node[elem[i].node[0]].y, node[elem[i].node[0]].z,
                          node[elem[i].node[2]].x, node[elem[i].node[2]].y, node[elem[i].node[2]].z);
            }
        }                       
    }

#endif

    void Interface3D::Trim(const UniformMesh3D& m, const int ksurf, const int kfilt, 
                           const int sval, const double fval, const int sgn, 
                           const double thresh, const char adjust)
    {
        int *mark = new int[ncount];
        
        // mark all nodes as inside (false) or outside (true)
        
        for (int i=0; i<ncount; ++i) {
            double temp = m.Interp(node[i].x, node[i].y, node[i].z, kfilt);
            mark[i] = sgn*(temp-fval) < 0.;
        }
        
        // throw out all elements which are completely outside
        
        for (int i=ecount-1; i>=0; --i) 
            if (mark[elem[i].node[0]] && mark[elem[i].node[1]] 
                && mark[elem[i].node[2]]) {
                --ecount;
                elem[i] = elem[ecount];
            }
        
        // throw out all loose nodes
        
        for (int i=0; i<ecount; ++i) {
            mark[elem[i].node[0]] = mark[elem[i].node[0]] ? -1 : 0;
            mark[elem[i].node[1]] = mark[elem[i].node[1]] ? -1 : 0;
            mark[elem[i].node[2]] = mark[elem[i].node[2]] ? -1 : 0;
        }
        int *xlate = new int[ncount];
        int ncnt = 0;
        for (int i=0; i<ncount; ++i) 
            switch(mark[i]) {
            case -1:
                mark[ncnt] = 1;
                node[ncnt] = node[i];
                xlate[i] = ncnt++;
                break;
            case 0:
                mark[ncnt] = 0;
                node[ncnt] = node[i];
                xlate[i] = ncnt++;
                break;
            case 1:
                break;
            };
        ncount = ncnt;
        for (int i=0; i<ecount; ++i) 
            for (int j=0; j<3; ++j)
                elem[i].node[j] = xlate[elem[i].node[j]];
        
        // push outside nodes back inside

        double grad[3], alpha;
        for (int i=0; i<ncount; ++i)
            if (mark[i]) {
                while(fabs(m.Interp(node[i].x,node[i].y,node[i].z,kfilt)-fval)>thresh
                      || fabs(m.Interp(node[i].x,node[i].y,node[i].z,ksurf)-sval)>thresh) {

                    // project onto filter surface

                    m.InterpGrad(node[i].x, node[i].y, node[i].z, kfilt, 
                                 grad[0], grad[1], grad[2]);
                    alpha = (fval-m.Interp(node[i].x,node[i].y,node[i].z,kfilt))
                        /(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
                    node[i].x += alpha*grad[0];
                    node[i].y += alpha*grad[1];
                    node[i].z += alpha*grad[2];
                                
                    // project onto interface surface
                                
                    m.InterpGrad(node[i].x, node[i].y, node[i].z, ksurf, 
                                 grad[0], grad[1], grad[2]);
                    alpha = (sval-m.Interp(node[i].x,node[i].y,node[i].z,ksurf))
                        /(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
                    node[i].x += alpha*grad[0];
                    node[i].y += alpha*grad[1];
                    node[i].z += alpha*grad[2];
                }
            }
        
        if (adjust) {
//#define DEBUG_I3D
#ifdef DEBUG_I3D
            PlotWindow3D win(true);
            win.SetRange(-2,2,-2,2,-2,2);
            win.SetView(1,0,100);
            Plot(win);
#endif
            int** web = new int*[ncount];
            for (int i=0; i<ncount; ++i) {
                web[i] = new int[20];
                web[i][0] = 0;
            }
            int ni, nj, k;
            for (int e=0; e<ecount; ++e) {
                for (int i=0; i<3; ++i) {
                    ni = elem[e].node[i];
                    for (int j=i+1; j<3; ++j) {
                        nj = elem[e].node[j];
                        if (mark[ni]) {
                            if (mark[nj]) {
                                for (k=1; k<=web[ni][0] && web[ni][k]!=nj; ++k) ;
                                if (k > web[ni][0]) {
                                    ++web[ni][0];
                                    web[ni][web[ni][0]] = nj;
                                    ++web[nj][0];
                                    web[nj][web[nj][0]] = ni;
                                }
                            } else {
                                for (k=1; k<=web[nj][0] && web[nj][k]!=ni; ++k) ;
                                if (k > web[nj][0]) {
                                    ++web[nj][0];
                                    web[nj][web[nj][0]] = ni;
                                }
                            }
                        } else {
                            if (mark[nj]) {
                                for (k=1; k<=web[ni][0] && web[ni][k]!=nj; ++k) ;
                                if (k > web[ni][0]) {
                                    ++web[ni][0];
                                    web[ni][web[ni][0]] = nj;
                                }
                            } else {
                                for (k=1; k<=web[ni][0] && web[ni][k]!=nj; ++k) ;
                                if (k > web[ni][0]) {
                                    ++web[ni][0];
                                    web[ni][web[ni][0]] = nj;
                                    ++web[nj][0];
                                    web[nj][web[nj][0]] = ni;
                                }
                            }
                        }
                    }
                }
            }
                
            double maxmove = 10000.;
            double nx, ny, nz;
            while (maxmove > thresh) {
                maxmove = 0.;
                for (int i=0; i<ncount; ++i) {
                    nx = 0.;
                    ny = 0.;
                    nz = 0.;
                    for (int j=1; j<=web[i][0]; ++j) {
                        nx += node[web[i][j]].x;
                        ny += node[web[i][j]].y;
                        nz += node[web[i][j]].z;
                    }
                    nx /= web[i][0];
                    ny /= web[i][0];
                    nz /= web[i][0];
                    if (mark[i]) {
                        while(fabs(m.Interp(nx,ny,nz,kfilt)-fval)>thresh
                              || fabs(m.Interp(nx,ny,nz,ksurf)-sval)>thresh) {

                            // project onto filter surface

                            m.InterpGrad(nx, ny, nz, kfilt, grad[0], grad[1], grad[2]);
                            alpha = (fval-m.Interp(nx,ny,nz,kfilt))
                                /(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
                            nx += alpha*grad[0];
                            ny += alpha*grad[1];
                            nz += alpha*grad[2];
                                                
                            // project onto interface surface
                                                
                            m.InterpGrad(nx, ny, nz, ksurf, grad[0], grad[1], grad[2]);
                            alpha = (sval-m.Interp(nx,ny,nz,ksurf))
                                /(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
                            nx += alpha*grad[0];
                            ny += alpha*grad[1];
                            nz += alpha*grad[2];
                        }
                    } else {
                        while(fabs(m.Interp(nx,ny,nz,ksurf)-sval)>thresh) {

                            // project onto interface surface
                                                
                            m.InterpGrad(nx, ny, nz, ksurf, grad[0], grad[1], grad[2]);
                            alpha = (sval-m.Interp(nx,ny,nz,ksurf))
                                /(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
                            nx += alpha*grad[0];
                            ny += alpha*grad[1];
                            nz += alpha*grad[2];
                        }
                    }
                    double move = sqrt(sqr(nx-node[i].x)+sqr(ny-node[i].y)+sqr(nz-node[i].z));
                    maxmove = maxmove > move ? maxmove : move;
                    node[i].x = nx; node[i].y = ny; node[i].z = nz;
                }
#ifdef DEBUG_I3D
                win.Clear();
                for (int e=0; e<ecount; ++e) {
                    win.Line(node[elem[e].node[0]].x,node[elem[e].node[0]].y,
                             node[elem[e].node[0]].z,node[elem[e].node[1]].x,
                             node[elem[e].node[1]].y,node[elem[e].node[1]].z);
                    win.Line(node[elem[e].node[2]].x,node[elem[e].node[2]].y,
                             node[elem[e].node[2]].z,node[elem[e].node[1]].x,
                             node[elem[e].node[1]].y,node[elem[e].node[1]].z);
                    win.Line(node[elem[e].node[0]].x,node[elem[e].node[0]].y,
                             node[elem[e].node[0]].z,node[elem[e].node[2]].x,
                             node[elem[e].node[2]].y,node[elem[e].node[2]].z);
                }
#endif
            }
                
            for (int i=0; i<ncount; ++i) delete[] web[i];
            delete[] web;
        }
        
        delete[] xlate;
        delete[] mark;
    }

    void Interface3D::Dump(const char* filename) const
    {
        std::ofstream s(filename);
        
        s << "Nodes:\n";
        for (int i=0; i<ncount; ++i) 
            s << i << ": (" << node[i].x << ',' << node[i].y << ',' << node[i].z << ")\n";
        s << "\nElements:\n";
        for (int i=0; i<ecount; ++i)
            s << i << ": (" << elem[i].node[0] << ',' << elem[i].node[1] << ',' 
              << elem[i].node[2] << ") = [(" << node[elem[i].node[0]].x << ','
              << node[elem[i].node[0]].y << ',' << node[elem[i].node[0]].z << "), ("
              << node[elem[i].node[1]].x << ',' << node[elem[i].node[1]].y << ','
              << node[elem[i].node[1]].z << "), (" << node[elem[i].node[2]].x << ','
              << node[elem[i].node[2]].y << ',' << node[elem[i].node[2]].z << ")]\n";
    }

}
