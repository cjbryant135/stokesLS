/*************************************************
 bicubicgrid.cpp
 
 Revision 1.0  2017/11/28  11:05:46  11:05:46  chopp (David Chopp)
 Initial revision
 
*************************************************/

#include "bicubicgrid.h"
#define DEBUG_BICUBICGRID
#ifdef DEBUG_BICUBICGRID
#include "plotwindow2d.h"
#endif

namespace levelset {

    BicubicGrid::BicubicGrid(const UniformMesh2D* themesh, const int k, const bool xperiodic, const bool yperiodic) 
        : mesh(themesh), grid(NULL), maxi(themesh->maxi), maxj(themesh->maxj), dx(themesh->dx), dy(themesh->dy), kval(k)
    {
        periodic[0] = xperiodic;
        periodic[1] = yperiodic;
        grid = new Bicubic*[maxi*maxj];
        for (int i=0; i<maxi*maxj; ++i)
            grid[i] = NULL;
        zero[0] = themesh->zero[0];
        zero[1] = themesh->zero[1];
    }

    BicubicGrid::~BicubicGrid(void)
    {
        if (grid) {
            for (int i=0; i<maxi*maxj; ++i)
                if (grid[i]) delete grid[i];
            delete[] grid;
        }
    }

    Bicubic* BicubicGrid::GetBicubic(const int ii, const int jj)
    {
        int i, j;
        if (periodic[0]) i=(ii+maxi)%maxi;
        else i = ii;
        if (periodic[1]) j=(jj+maxj)%maxj;
        else j = jj;
        Bicubic* p = grid[Index(i,j)];
        if (p == NULL) {
            p = grid[Index(i,j)] = new Bicubic;
            double data[4][4];
            if (periodic[0]) {
                if (periodic[1]) {
                    for (int k=0; k<4; ++k)
                        for (int l=0; l<4; ++l)
                            data[k][l] = mesh->data_((i+k-1+maxi)%maxi, (j+l-1+maxj)%maxj, kval);
                } else {
                    for (int k=0; k<4; ++k)
                        for (int l=0; l<4; ++l)
                            data[k][l] = mesh->data_((i+k-1+maxi)%maxi, j+l-1, kval);
                }
            } else {
                if (periodic[1]) {
                    for (int k=0; k<4; ++k)
                        for (int l=0; l<4; ++l)
                            data[k][l] = mesh->data_(i+k-1, (j+l-1+maxj)%maxj, kval);
                } else {
                    for (int k=0; k<4; ++k)
                        for (int l=0; l<4; ++l)
                            data[k][l] = mesh->data_(i+k-1, j+l-1, kval);
                }
            }

            p->BuildwDeriv(data[0][0], data[1][0], data[2][0], data[3][0],
                           data[0][1], data[1][1], data[2][1], data[3][1],
                           data[0][2], data[1][2], data[2][2], data[3][2],
                           data[0][3], data[1][3], data[2][3], data[3][3],
                           dx, dy);
        }
        return p;
    }

    Bicubic* BicubicGrid::GetBicubic(const double x, const double y)
    {
        return GetBicubic(I(x), J(y));
    }

    double BicubicGrid::operator()(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->F(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dx(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dx(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dy(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dy(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dxx(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dxx(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dxy(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dxy(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dyy(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dyy(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dxxx(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dxxx(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dxxy(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dxxy(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dxyy(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dxyy(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::Dyyy(const double x, const double y) 
    {
        Bicubic* p = GetBicubic(x,y);
        return p->Dyyy(x-I(x)*dx,y-J(y)*dy);
    }
    
    double BicubicGrid::LocalDist(const double x, const double y,
                                  double& ax, double& ay, char& clean)
    {
        ax = x;
        ay = y;
        double tol = 1.e-3;
        double p, px, py, tax, tay;
        double delta[4] = {1000.,1000.,1000.,1000.};
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
                
        while (resid > tol*dx*dy && (itcount < 20 || resid < oldresid)) {
                        
            // Find zero
                        
            p = F(ax,ay);
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            delta[0] = -p*px/(px*px+py*py);
            delta[1] = -p*py/(px*px+py*py);
            tax = ax+delta[0];
            tay = ay+delta[1];
                        
            // Find min distance
                        
            p = F(tax,tay);
            px = Dx(tax,tay);
            py = Dy(tax,tay);
            double denom = px*px+py*py;
            delta[2] = py*(px*(tay-y)-py*(tax-x))/denom;
            delta[3] = -px*(px*(tay-y)-py*(tax-x))/denom;
                        
            ax += delta[0]+delta[2];
            ay += delta[1]+delta[3];
            ++itcount;
            oldresid = resid;
            resid = sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3]));
        }
                
        //    clean = itcount < 20;
        clean = resid < tol*dx*dy && ax >= 0. && ax <= dx && ay >= 0. && ay <= dy;
        
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }

#if 0
    double BicubicGrid::LocalDistNewton(const double x, const double y, 
                                        double& ax, double& ay, char& clean)
    {
        double tol = 1.e-3;
        double p, px, py, pxx, pxy, pyy;
#if 0
        px = Dx(x,y);
        py = Dy(x,y);
        ax = x + sign(px)*dx/2.;
        ay = y + sign(py)*dy/2.;
#endif
        double Jinv[2][2];
        double det;
        double delta[2];
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
#ifdef DEBUG_BICUBICGRID
        bool usewindow = false;
        int winwidth=1;
        PlotWindow2D window(usewindow, usewindow);
        if (usewindow) {
            window.SetRange(x-winwidth*dx,x+winwidth*dx,y-winwidth*dy,y+winwidth*dy);
            for (int i=-winwidth+1; i<=winwidth-1; ++i) {
                window.Line(x-winwidth*dx,y+i*dy,x+winwidth*dx,y+i*dy);
                window.Line(x+i*dx,y-winwidth*dy,x+i*dx,y+winwidth*dy);
            }
            double ldx = dx/50;
            double ldy = dy/50;
            double val[2][2];
            for (int i=0; i<winwidth*100; ++i) {
                double lx = x-winwidth*dx+i*ldx;
                for (int j=0; j<winwidth*100; ++j) {
                    double ly = y-winwidth*dy+j*ldy;
                    for (int ii=0; ii<=1; ++ii)
                        for (int jj=0; jj<=1; ++jj)
                            val[ii][jj] = F(lx+ii*ldx,ly+jj*ldy);
                    int score = (val[0][0] > 0. ? 1 : 0) + (val[1][0] > 0. ? 2 : 0)
                        + (val[0][1] > 0. ? 4 : 0) + (val[1][1] > 0. ? 8 : 0);
                    switch (score) {
                    case 1:
                    case 14:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]));
                        break;
                    case 2:
                    case 13:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    case 3:
                    case 12:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    case 4:
                    case 11:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 5:
                    case 10:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 6:
                    case 9:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 7:
                    case 8:
                        window.Line(lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy,
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    }
                }
            }
            window.Dot(ax,ay);
        }
#endif

        bool doublecheck = false;
        while (resid > tol*dx*dy && (itcount < 20 || resid < oldresid)) {
                        
            // Find zero
                        
            p = F(ax,ay);
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            pxx = Dxx(ax,ay);
            pxy = Dxy(ax,ay);
            pyy = Dyy(ax,ay);
            det = px*(pxy*(ay-y)+px-pyy*(ax-x)) + py*(pxy*(ax-x)+py-pxx*(ay-y));
            Jinv[0][0] = (pxy*(ay-y)+px-pyy*(ax-x))/det;
            Jinv[0][1] = -py/det;
            Jinv[1][0] = (pxy*(ax-x)+py-pxx*(ay-y))/det;
            Jinv[1][1] = px/det;
            delta[0] = -Jinv[0][0]*p-Jinv[0][1]*(px*(ay-y)-py*(ax-x));
            delta[1] = -Jinv[1][0]*p-Jinv[1][1]*(px*(ay-y)-py*(ax-x));
            if (fabs(delta[0]) > dx/2.) {
                delta[0] *= dx/2./fabs(delta[0]);
                delta[1] *= dx/2./fabs(delta[0]);
            }
            if (fabs(delta[1]) > dy/2.) {
                delta[0] *= dy/2./fabs(delta[1]);
                delta[1] *= dy/2./fabs(delta[1]);
            }
                                                
            ax += delta[0];
            ay += delta[1];
#if 0
            if (ax < 0.) ax = ax+maxi*dx;
            if (ax > X(maxi-1)) ax = ax-maxi*dx;
            if (ay < 0.) ay = ay+maxj*dy;
            if (ay > Y(maxj-1)) ay = ay-maxj*dy;
#endif
#ifdef DEBUG_BICUBICGRID
            if (usewindow)
                window.Dot(ax,ay);
#endif
            ++itcount;
            oldresid = resid;
            resid = sqrt(sqr(delta[0])+sqr(delta[1]));
            if (!doublecheck && itcount >= 20 && resid > tol*dx*dy && oldresid < resid) {
                doublecheck = true;
                ax = x;
                ay = y;
                itcount = 0;
            }
        }
                
        //    clean = itcount < 20;
        clean = resid < tol*dx*dy;
        
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }
#else

    // New version using Lagrange multiplier method
    double BicubicGrid::LocalDistNewton(const double x, const double y, 
                                        double& ax, double& ay, char& clean)
    {
        double tol = 1.e-3;
        double p, px, py, pxx, pxy, pyy;
        double lambda = 0.;
        double Jinv[3][3];
        double det;
        double delta[3];
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
#ifdef DEBUG_BICUBICGRID
        bool usewindow = false;
        int winwidth=2;
        PlotWindow2D window(usewindow, false);
        if (usewindow) {
            window.SetRange(x-winwidth*dx,x+winwidth*dx,y-winwidth*dy,y+winwidth*dy);
            for (int i=-winwidth+1; i<=winwidth-1; ++i) {
                window.Line(x-winwidth*dx,y+i*dy,x+winwidth*dx,y+i*dy);
                window.Line(x+i*dx,y-winwidth*dy,x+i*dx,y+winwidth*dy);
            }
            double ldx = dx/50;
            double ldy = dy/50;
            double val[2][2];
            for (int i=0; i<winwidth*100; ++i) {
                double lx = x-winwidth*dx+i*ldx;
                for (int j=0; j<winwidth*100; ++j) {
                    double ly = y-winwidth*dy+j*ldy;
                    for (int ii=0; ii<=1; ++ii)
                        for (int jj=0; jj<=1; ++jj)
                            val[ii][jj] = F(lx+ii*ldx,ly+jj*ldy);
                    int score = (val[0][0] > 0. ? 1 : 0) + (val[1][0] > 0. ? 2 : 0)
                        + (val[0][1] > 0. ? 4 : 0) + (val[1][1] > 0. ? 8 : 0);
                    switch (score) {
                    case 1:
                    case 14:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]));
                        break;
                    case 2:
                    case 13:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    case 3:
                    case 12:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    case 4:
                    case 11:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 5:
                    case 10:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 6:
                    case 9:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 7:
                    case 8:
                        window.Line(lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy,
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    }
                }
            }
            window.Dot(ax,ay);
        }
#endif

        bool doublecheck = false;
        while (resid > tol*dx*dy && (itcount < 20 || resid < oldresid)) {
                        
            // Find zero
                        
            p = F(ax,ay);
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            pxx = lambda*Dxx(ax,ay);
            pxy = lambda*Dxy(ax,ay);
            pyy = lambda*Dyy(ax,ay);
            det = -px*px*pyy + 2*px*py*pxy - py*py*pxx - 2*px*px - 2*py*py;
            Jinv[0][0] = -py*py/det;
            Jinv[0][1] = px*py/det;
            Jinv[0][2] = (-px*pyy + pxy*py - 2*px)/det;
            Jinv[1][0] = Jinv[0][1];
            Jinv[1][1] = -px*px/det;
            Jinv[1][2] = (px*pxy - pxx*py - 2*py)/det;
            Jinv[2][0] = Jinv[0][2];
            Jinv[2][1] = Jinv[1][2];
            Jinv[2][2] = (pxx*pyy - pxy*pxy + 2*pxx + 2*pyy + 4)/det;
            delta[0] = -Jinv[0][0]*(2*(ax-x)+lambda*px)-Jinv[0][1]*(2*(ay-y)+lambda*py)-Jinv[0][2]*p;
            delta[1] = -Jinv[1][0]*(2*(ax-x)+lambda*px)-Jinv[1][1]*(2*(ay-y)+lambda*py)-Jinv[1][2]*p;
            delta[2] = -Jinv[2][0]*(2*(ax-x)+lambda*px)-Jinv[2][1]*(2*(ay-y)+lambda*py)-Jinv[2][2]*p;
#if 1
            if (fabs(delta[0]) > dx/2.) {
                delta[0] *= dx/2./fabs(delta[0]);
                delta[1] *= dx/2./fabs(delta[0]);
                delta[2] *= dx/2./fabs(delta[0]);
            }
            if (fabs(delta[1]) > dy/2.) {
                delta[0] *= dy/2./fabs(delta[1]);
                delta[1] *= dy/2./fabs(delta[1]);
                delta[2] *= dy/2./fabs(delta[1]);
            }
#endif                                          
            ax += delta[0];
            ay += delta[1];
            lambda += delta[2];
#ifdef DEBUG_BICUBICGRID
            if (usewindow)
                window.Dot(ax,ay);
#endif
            ++itcount;
            oldresid = resid;
            resid = sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2]));
            if (!doublecheck && itcount >= 20 && resid > tol*dx*dy && oldresid < resid) {
                doublecheck = true;
                ax = x;
                ay = y;
                itcount = 0;
            }
        }
                
        //    clean = itcount < 20;
        clean = resid < tol*dx*dy;
        
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }
#endif

    double DistToSeg(const double x0, const double y0, const double x1, const double y1, 
                     const double x, const double y, double& ax, double& ay)
    {
        double dx[2], dy[2];
        dx[0] = x1-x0;
        dx[1] = x-x0;
        dy[0] = y1-y0;
        dy[1] = y-y0;
        double t = (dx[0]*dx[1]+dy[0]*dy[1])/(sqr(dx[0])+sqr(dy[0]));
        if (t < 0.) {
            ax = x0; ay = y0; 
            return sqrt(sqr(dx[1])+sqr(dy[1]));
        } else if (t > 1.) {
            ax = x1; ay = y1; 
            return sqrt(sqr(x-x1)+sqr(y-y1));
        } else {
            double t1 = dx[1]*dx[0]+dy[1]*dy[0];
            double t2 = (x1-x)*dx[0]+(y1-y)*dy[0];
            double denom = sqr(dx[0])+sqr(dy[0]);
            ax = (t1*x1+t2*x0)/denom;
            ay = (t1*y1+t2*y0)/denom;
            return sqrt(sqr(ax-x)+sqr(ay-y));
        }
    }

    double BoxDist(const double f00, const double f10, const double f01, 
                   const double f11, const double dx, const double dy,
                   double& ax, double& ay)
    {
        double dist;
        int score = (f00 >= 0. ? 1 : 0) + (f10 >= 0. ? 2 : 0) 
            + (f01 >= 0. ? 4 : 0) + (f11 >= 0. ? 8 : 0);
        if (score > 7) score = 15-score;
        switch(score) {
        case 0:
            dist = DBL_MAX; 
            ax = dx; 
            ay = dy;
            break;
        case 1:
        case 6:
            dist = DistToSeg(-dx*f00/(f10-f00),0.,0.,-dy*f00/(f01-f00),0.,0.,ax,ay);
            break;
        case 2:
            ax = -dx*f00/(f10-f00);
            ay = 0.;
            dist = ax;
            break;
        case 3:
            dist = DistToSeg(0.,-dy*f00/(f01-f00),dx,-dy*f10/(f11-f10),0.,0.,ax,ay);
            break;
        case 4:
            ax = 0.;
            ay = -dy*f00/(f01-f00);
            dist = ay;
            break;
        case 5:
            dist = DistToSeg(-dx*f00/(f10-f00),0.,-dx*f01/(f11-f01),dy,0.,0.,ax,ay);
            break;
        case 7:
            dist = DistToSeg(dx,-dy*f10/(f11-f10),-dx*f01/(f11-f01),dy,0.,0.,ax,ay);
            break;
        }
    }
      


    double BicubicGrid::LocalDistDirect(const double x, const double y,
                                        double& ax, double& ay)
    {
        // check all four surrounding boxes
        int sgn[3][3];
        double dist = DBL_MAX;
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                sgn[i][j] = (F(x+(i-1)*dx,y+(j-1)*dy) >= 0. ? 1 : 0);
        // check lower left
        int score = sgn[0][0]+sgn[1][0]+sgn[0][1]+sgn[1][1];
        if (score > 0 && score < 4) {
            double tax, tay;
            double tempdist = BoxDist(F(x,y),F(x-dx,y),F(x,y-dy),F(x-dx,y-dy),dx,dy,tax,tay);
            if (tempdist < dist) {
                dist = tempdist;
                ax = x-tax;
                ay = y-tay;
            }
        }
        // check lower right
        score = sgn[1][0]+sgn[2][0]+sgn[1][1]+sgn[2][1];
        if (score > 0 && score < 4) {
            double tax, tay;
            double tempdist = BoxDist(F(x,y),F(x+dx,y),F(x,y-dy),F(x+dx,y-dy),dx,dy,tax,tay);
            if (tempdist < dist) {
                dist = tempdist;
                ax = x+tax;
                ay = y-tay;
            }
        }
        // check upper left
        score = sgn[0][1]+sgn[1][1]+sgn[0][2]+sgn[1][2];
        if (score > 0 && score < 4) {
            double tax, tay;
            double tempdist = BoxDist(F(x,y),F(x-dx,y),F(x,y+dy),F(x-dx,y+dy),dx,dy,tax,tay);
            if (tempdist < dist) {
                dist = tempdist;
                ax = x-tax;
                ay = y+tay;
            }
        }
        // check upper right
        score = sgn[1][1]+sgn[2][1]+sgn[1][2]+sgn[2][2];
        if (score > 0 && score < 4) {
            double tax, tay;
            double tempdist = BoxDist(F(x,y),F(x+dx,y),F(x,y+dy),F(x+dx,y+dy),dx,dy,tax,tay);
            if (tempdist < dist) {
                dist = tempdist;
                ax = x+tax;
                ay = y+tay;
            }
        }
        return dist;
    }


    double BicubicGrid::LocalDistBisection(const double x, const double y, 
                                           double& ax, double& ay, char& clean)
    {
        double tol = 1.e-3;
        double p, px, py, pxx, pxy, pyy;
        double Jinv[2][2];
        double det;
        double delta[2];
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
#ifdef DEBUG_BICUBICGRID
        bool usewindow = false;
        int winwidth=1;
        PlotWindow2D window(usewindow, usewindow);
        if (usewindow) {
            window.SetRange(x-winwidth*dx,x+winwidth*dx,y-winwidth*dy,y+winwidth*dy);
            for (int i=-winwidth+1; i<=winwidth-1; ++i) {
                window.Line(x-winwidth*dx,y+i*dy,x+winwidth*dx,y+i*dy);
                window.Line(x+i*dx,y-winwidth*dy,x+i*dx,y+winwidth*dy);
            }
            double ldx = dx/50;
            double ldy = dy/50;
            double val[2][2];
            for (int i=0; i<winwidth*100; ++i) {
                double lx = x-winwidth*dx+i*ldx;
                for (int j=0; j<winwidth*100; ++j) {
                    double ly = y-winwidth*dy+j*ldy;
                    for (int ii=0; ii<=1; ++ii)
                        for (int jj=0; jj<=1; ++jj)
                            val[ii][jj] = F(lx+ii*ldx,ly+jj*ldy);
                    int score = (val[0][0] > 0. ? 1 : 0) + (val[1][0] > 0. ? 2 : 0)
                        + (val[0][1] > 0. ? 4 : 0) + (val[1][1] > 0. ? 8 : 0);
                    switch (score) {
                    case 1:
                    case 14:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]));
                        break;
                    case 2:
                    case 13:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    case 3:
                    case 12:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    case 4:
                    case 11:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 5:
                    case 10:
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 6:
                    case 9:
                        window.Line(lx, ly-val[0][0]*ldy/(val[0][1]-val[0][0]),
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        window.Line(lx-val[0][0]*ldx/(val[1][0]-val[0][0]), ly,
                                    lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy);
                        break;
                    case 7:
                    case 8:
                        window.Line(lx-val[0][1]*ldx/(val[1][1]-val[0][1]), ly+ldy,
                                    lx+ldx, ly-val[1][0]*ldy/(val[1][1]-val[1][0]));
                        break;
                    }
                }
            }
            window.Dot(ax,ay);
        }
#endif

        // collapse to interface
        p = F(ax,ay);
        while (fabs(p) > tol*dx*dy) {
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            ax -= p*px/(px*px+py*py);
            ay -= p*py/(px*px+py*py);
            p = F(ax,ay);
#ifdef DEBUG_BICUBICGRID
            window.Dot(ax,ay);
#endif
        }
        px = Dx(ax,ay);
        py = Dy(ax,ay);
        double cp = (ax-x)*py-(ay-y)*px;
        double ds = min(dx,dy)/100.;
        double ax2 = ax-ds*py/sqrt(px*px+py*py);
        double ay2 = ay+ds*px/sqrt(px*px+py*py);
        double p2 = F(ax2,ay2);
        double px2, py2;
        while (fabs(p2) > tol*dx*dy) {
            px2 = Dx(ax2,ay2);
            py2 = Dy(ax2,ay2);
            ax2 -= p2*px2/(px2*px2+py2*py2);
            ay2 -= p2*py2/(px2*px2+py2*py2);
            p2 = F(ax2,ay2);
#ifdef DEBUG_BICUBICGRID
            window.Dot(ax2,ay2);
#endif
        }
        px2 = Dx(ax2,ay2);
        py2 = Dy(ax2,ay2);
        double cp2 = (ax2-x)*py2-(ay2-y)*px2;
        while (cp2*cp > 0) {
            ax2 = ax2 + ds*py2/sqrt(px2*px2+py2*py2);
            ay2 = ay2 - ds*px2/sqrt(px2*px2+py2*py2);
#ifdef DEBUG_BICUBICGRID
            window.Dot(ax2,ay2);
#endif
            p2 = F(ax2,ay2);
            while (fabs(p2) > tol*dx) {
                px2 = Dx(ax2,ay2);
                py2 = Dy(ax2,ay2);
                ax2 -= p2*px2/(px2*px2+py2*py2);
                ay2 -= p2*py2/(px2*px2+py2*py2);
                p2 = F(ax2,ay2);
#ifdef DEBUG_BICUBICGRID
                window.Dot(ax2,ay2);
#endif
            }
            px2 = Dx(ax2,ay2);
            py2 = Dy(ax2,ay2);
            cp2 = (ax2-x)*py2-(ay2-y)*px2;
        }

        // now do bisection

        while (fabs(ax-ax2)+fabs(ay-ay2) > tol*dx*dy) {
            double ax3 = (ax+ax2)/2.;
            double ay3 = (ay+ay2)/2.;
#ifdef DEBUG_BICUBICGRID
            window.Dot(ax3,ay3);
#endif
            double p3 = F(ax3,ay3);
            double px3;
            double py3;
            while (fabs(p3) > tol*dx*dy) {
                px3 = Dx(ax3,ay3);
                py3 = Dy(ax3,ay3);
                ax3 -= p3*px3/(px3*px3+py3*py3);
                ay3 -= p3*py3/(px3*px3+py3*py3);
                p3 = F(ax3,ay3);
#ifdef DEBUG_BICUBICGRID
                window.Dot(ax3,ay3);
#endif
            }
            px3 = Dx(ax3,ay3);
            py3 = Dy(ax3,ay3);
            double cp3 = (ax3-x)*py3-(ay3-y)*px3;
            if (sign(cp3) == sign(cp)) {
                ax = ax3;
                ay = ay3;
                cp = cp3;
            } else {
                ax2 = ax3;
                ay2 = ay3;
                cp2 = cp3;
            }
        }
        //    clean = resid < tol*dx*dy;
        clean = 1;
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }

    void BicubicGrid::FollowToBdry(const int dir, const double x, const double y, double& ax, double &ay)
    {
        Bicubic* p = GetBicubic(x,y);
        int i = I(x);
        int j = J(y);
        double lx = x-X(i);
        double ly = y-Y(j);
        p->FollowToBdry(dir, lx, ly, ax, ay);
    }
}
  
  
