/*************************************************
 bicubic.cpp
 
 $Header: bicubic.cpp,v 1.1 2000/05/31 11:05:46 chopp Exp $
 
 $Log:  bicubic.cpp,v $
 Revision 1.1  2000/05/31  11:05:46  11:05:46  chopp (David Chopp)
 Initial revision
 
*************************************************/

#include "bicubic.h"
//#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stack>
#include "lapack.h"

using namespace std;

namespace levelset {
        
    Bicubic::Bicubic(const Bicubic& p)
    {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                a[i][j] = p.a[i][j];
        wx = p.wx;
        wy = p.wy;
    }
        
        
    void Bicubic::Build(const double f11, const double f21, const double f12,
                        const double f22, const double fx11, const double fx21,
                        const double fx12, const double fx22, const double fy11,
                        const double fy21, const double fy12, const double fy22,
                        const double fxy11, const double fxy21, const double fxy12,
                        const double fxy22, const double dx, const double dy)
    {
        double f[16];
        f[0] = f11;
        f[1] = f21;
        f[2] = f12;
        f[3] = f22;
        f[4] = fx11*dx;
        f[5] = fx21*dx;
        f[6] = fx12*dx;
        f[7] = fx22*dx;
        f[8] = fy11*dy;
        f[9] = fy21*dy;
        f[10] = fy12*dy;
        f[11] = fy22*dy;
        f[12] = fxy11*dx*dy;
        f[13] = fxy21*dx*dy;
        f[14] = fxy12*dx*dy;
        f[15] = fxy22*dx*dy;
                
        a[0][0] = f[0];
        a[0][1] = f[8]/dy;
        a[0][2] = (-3*f[0]+3*f[2]-2*f[8]-f[10])/dy/dy;
        a[0][3] = (2*f[0]-2*f[2]+f[8]+f[10])/dy/dy/dy;
        a[1][0] = f[4]/dx;                               
        a[1][1] = f[12]/dx/dy;         
        a[1][2] = (-3*f[4]+3*f[6]-2*f[12]-f[14])/dx/dy/dy; 
        a[1][3] = (2*f[4]-2*f[6]+f[12]+f[14])/dx/dy/dy/dy;
        a[2][0] = (-3*f[0]+3*f[1]-2*f[4]-f[5])/dx/dx;                             
        a[2][1] = (-3*f[8]+3*f[9]-2*f[12]-f[13])/dx/dx/dy;   
        a[2][2] =  (9*f[0]-9*f[1]-9*f[2]+9*f[3]+6*f[4]+3*f[5]-6*f[6]-3*f[7]+6*f[8]
                    -6*f[9]+3*f[10]-3*f[11]+4*f[12]+2*f[13]+2*f[14]+f[15])
            /dx/dx/dy/dy;
        a[2][3] = (-6*f[0]+6*f[1]+6*f[2]-6*f[3]-4*f[4]-2*f[5]+4*f[6]+2*f[7]-3*f[8] 
                   +3*f[9]-3*f[10]+3*f[11]-2*f[12]-f[13]-2*f[14]-f[15])
            /dx/dx/dy/dy/dy;
        a[3][0] =  (2*f[0]-2*f[1]+f[4]+f[5])/dx/dx/dx;
        a[3][1] =  (2*f[8]-2*f[9]+f[12]+f[13])/dx/dx/dx/dy;
        a[3][2] = (-6*f[0]+6*f[1]+6*f[2]-6*f[3]-3*f[4]-3*f[5]+3*f[6]+3*f[7]-4*f[8]
                   +4*f[9]-2*f[10]+2*f[11]-2*f[12]-2*f[13]-f[14]-f[15])
            /dx/dx/dx/dy/dy;
        a[3][3] =  (4*f[0]-4*f[1]-4*f[2]+4*f[3]+2*f[4]+2*f[5]-2*f[6]-2*f[7]+2*f[8]
                    -2*f[9]+2*f[10]-2*f[11]+f[12]+f[13]+f[14]+f[15])
            /dx/dx/dx/dy/dy/dy;
        wx = dx;
        wy = dy;
    }
        
    void Bicubic::BuildwDeriv(const double f00, const double f10, const double f20,
                              const double f30, const double f01, const double f11,
                              const double f21, const double f31, const double f02,
                              const double f12, const double f22, const double f32,
                              const double f03, const double f13, const double f23,
                              const double f33, const double dx, const double dy)
    {
#if 0
        Build(f11, f21, f12, f22, (f21-f01)/2./dx, (f31-f11)/2./dx, (f22-f02)/2./dx,
              (f32-f12)/2./dx, (f12-f10)/2./dy, (f22-f20)/2./dy, (f13-f11)/2./dy,
              (f23-f21)/2./dy, (f22-f20-f02+f00)/4./dx/dy, (f32-f30-f12+f10)/4./dx/dy,
              (f23-f21-f03+f01)/4./dx/dy, (f33-f31-f13+f11)/4./dx/dy, dx, dy);
#else
        const double A[16][16] = {{  0.,   0.,   0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
                                  {  0.,   0.,   0.,   0., -0.5,   0.,  0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
                                  {  0.,   0.,   0.,   0.,   1., -2.5,   2., -0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
                                  {  0.,   0.,   0.,   0., -0.5,  1.5, -1.5,  0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
                                  {  0., -0.5,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  0.5,   0.,   0.,   0.,   0.,   0.,   0.},
                                  {0.25,   0.,-0.25,   0.,   0.,   0.,   0.,   0.,-0.25,   0., 0.25,   0.,   0.,   0.,   0.,   0.},
                                  {-0.5, 1.25,  -1., 0.25,   0.,   0.,   0.,   0.,  0.5,-1.25,   1.,-0.25,   0.,   0.,   0.,   0.},
                                  {0.25,-0.75, 0.75,-0.25,   0.,   0.,   0.,   0.,-0.25, 0.75,-0.75, 0.25,   0.,   0.,   0.,   0.},
                                  {  0.,   1.,   0.,   0.,   0., -2.5,   0.,   0.,   0.,   2.,   0.,   0.,   0., -0.5,   0.,   0.},
                                  {-0.5,   0.,  0.5,   0., 1.25,   0.,-1.25,   0.,  -1.,   0.,   1.,   0., 0.25,   0.,-0.25,   0.},
                                  {  1., -2.5,   2., -0.5, -2.5, 6.25,  -5., 1.25,   2.,  -5.,   4.,  -1., -0.5, 1.25,  -1., 0.25},
                                  {-0.5,  1.5, -1.5,  0.5, 1.25,-3.75, 3.75,-1.25,  -1.,   3.,  -3.,   1., 0.25,-0.75, 0.75,-0.25},
                                  {  0., -0.5,   0.,   0.,   0.,  1.5,   0.,   0.,   0., -1.5,   0.,   0.,   0.,  0.5,   0.,   0.},
                                  {0.25,   0.,-0.25,   0.,-0.75,   0., 0.75,   0., 0.75,   0.,-0.75,   0.,-0.25,   0., 0.25,   0.},
                                  {-0.5, 1.25,  -1., 0.25,  1.5,-3.75,   3.,-0.75, -1.5, 3.75,  -3., 0.75,  0.5,-1.25,   1.,-0.25},
                                  {0.25,-0.75, 0.75,-0.25,-0.75, 2.25,-2.25, 0.75, 0.75,-2.25, 2.25,-0.75,-0.25, 0.75,-0.75, 0.25}};
                
        double f[16] = {f00, f01, f02, f03, f10, f11, f12, f13, f20, f21, f22, f23, f30, f31, f32, f33};
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j) {
                a[i][j] = 0.;
                for (int k=0; k<16; ++k)
                    a[i][j] += A[4*i+j][k]*f[k]/pow(dx,i)/pow(dy,j);
            }
        wx = dx;
        wy = dy;
#endif
    }
        
        
#ifdef PARTIAL_BICUBIC
    void Bicubic::BuildwDeriv(const double f00, const double f10, const double f20,
                              const double f30, const double f01, const double f11,
                              const double f21, const double f31, const double f02,
                              const double f12, const double f22, const double f32,
                              const double f03, const double f13, const double f23,
                              const double f33, const double dx, const double dy, 
                              const char* valid)
    {
        const double A[16][16] = {{1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1},
                                  {1,0,0,0,-1,0,0,0,1,0,0,0,1,0,0,0},
                                  {1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,1,1,1},
                                  {1,2,4,6,-1,-2,-4,-6,1,2,4,6,1,2,4,6},
                                  {1,-1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
                                  {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                  {1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
                                  {1,2,4,6,0,0,0,0,0,0,0,0,0,0,0,0},
                                  {1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1},
                                  {1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0},
                                  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
                                  {1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6},
                                  {1,-1,1,1,2,-2,2,2,4,-4,4,4,6,-6,6,6},
                                  {1,0,0,0,2,0,0,0,4,0,0,0,6,0,0,0},
                                  {1,1,1,1,2,2,2,2,4,4,4,4,6,6,6,6},
                                  {1,2,4,6,2,4,8,12,4,8,16,24,6,12,24,36}};
                
        double B[16];
        for (int i=0; i<16; ++i) B[i] = pow(dx,i%4)*pow(dy,i/4);
                
        double b[16];
        b[0] = f00; b[1] = f10; b[2] = f20; b[3] = f30;
        b[4] = f01; b[5] = f11; b[6] = f21; b[7] = f31;
        b[8] = f02; b[9] = f12; b[10] = f22; b[11] = f32;
        b[12] = f03; b[13] = f13; b[14] = f23; b[15] = f33;
                
        double M[16][16];       
        int vcount = 0;
        for (int r=0; r<16; ++r) {
            if (valid[r]) {
                ++vcount;
                for (int c=0; c<16; ++c)
                    M[r][c] = A[r][c]*B[c];
            } else {
                b[r] = 0.;
                for (int c=0; c<16; ++c)
                    M[r][c] = 0.;
            }
        }
                
        double U[16][16], V[16][16], S[16];
        double temp[128];
        long int info;
        dgesvd('A', 'A', 16, 16, &(M[0][0]), 16, S, &(U[0][0]), 16, 
               &(V[0][0]), 16, temp, 128, info);
                
        for (int i=0; i<16; ++i) {
            temp[i] = 0.;
            for (int j=0; j<16; ++j)
                temp[i] += V[j][i]*b[j];
        }
                
        for (int i=0; i<vcount; ++i)
            temp[i] /= S[i];
        for (int i=vcount; i<16; ++i)
            temp[i] = 0.;
                
        for (int i=0; i<16; ++i) {
            a[i%4][i/4] = 0.;
            for (int j=0; j<16; ++j) 
                a[i%4][i/4] += U[j][i]*temp[j];
        }
        wx = dx;
        wy = dy;
    }
#endif
        
        
    double Bicubic::operator()(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j)
                ans += a[i][j]*pow(x,i)*pow(y,j);
                
        return ans;
    }
        
    double Bicubic::Dx(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=1; i<4; ++i) 
            for (int j=0; j<4; ++j)
                ans += i*a[i][j]*pow(x,i-1)*pow(y,j);
                
        return ans;
    }
        
    double Bicubic::Dy(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=0; i<4; ++i) 
            for (int j=1; j<4; ++j)
                ans += j*a[i][j]*pow(x,i)*pow(y,j-1);
                
        return ans;
    }
        
    double Bicubic::Dxx(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=2; i<4; ++i) 
            for (int j=0; j<4; ++j)
                ans += i*(i-1)*a[i][j]*pow(x,i-2)*pow(y,j);
                
        return ans;
    }
        
    double Bicubic::Dxy(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=1; i<4; ++i) 
            for (int j=1; j<4; ++j)
                ans += i*j*a[i][j]*pow(x,i-1)*pow(y,j-1);
                
        return ans;
    }
        
    double Bicubic::Dyy(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=0; i<4; ++i) 
            for (int j=2; j<4; ++j)
                ans += j*(j-1)*a[i][j]*pow(x,i)*pow(y,j-2);
                
        return ans;
    }
        
    double Bicubic::Dxxx(const double x, const double y) const
    {
        double ans = 0;
                
        for (int j=0; j<4; ++j)
            ans += 6*a[3][j]*pow(y,j);
                
        return ans;
    }
        
    double Bicubic::Dxxy(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=2; i<4; ++i) 
            for (int j=1; j<4; ++j)
                ans += i*(i-1)*j*a[i][j]*pow(x,i-2)*pow(y,j-1);
                
        return ans;
    }
        
    double Bicubic::Dxyy(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=1; i<4; ++i) 
            for (int j=2; j<4; ++j)
                ans += i*j*(j-1)*a[i][j]*pow(x,i-1)*pow(y,j-2);
                
        return ans;
    }
        
    double Bicubic::Dyyy(const double x, const double y) const
    {
        double ans = 0;
                
        for (int i=0; i<4; ++i) 
            ans += 6*a[i][3]*pow(x,i);
                
        return ans;
    }
        
    Bicubic Bicubic::Dx(void)
    {
        Bicubic d;
                
        for (int j=0; j<4; ++j) {
            d.a[3][j] = 0.;
            for (int i=1; i<4; ++i)
                d.a[i-1][j] = i*a[i][j];
        }
                
        return d;
    }
        
    Bicubic Bicubic::Dy(void)
    {
        Bicubic d;
                
        for (int i=0; i<4; ++i) {
            d.a[i][3] = 0.;
            for (int j=1; j<4; ++j)
                d.a[i][j-1] = j*a[i][j];
        }
                
        return d;
    }
        
    inline double sqr(const double x) {return x*x;}

  class Box {
  public:
    double left, right, bottom, top;
    Box(double l=0., double r=0., double b=0., double t=0.) : left(l), right(r), bottom(b), top(t) {}
    Box Subbox(char isleft, char isbottom) {
      Box newbox(left, right, bottom, top);
      if (isleft)
	newbox.right = (left+right)/2.;
      else
	newbox.left = (left+right)/2.;
      if (isbottom)
	newbox.top = (bottom+top)/2.;
      else
	newbox.bottom = (bottom+top)/.2;
      return newbox;
    }
  };
        
    double Bicubic::LocalDist(const double x, const double y, 
                              double& ax, double& ay, char& clean) const
    {
        ax = x;
        ay = y;
        double tol = 1.e-3;
        double p, px, py, tax, tay;
        double delta[4] = {1000.,1000.,1000.,1000.};
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
                
        while (resid > tol*wx*wy && (itcount < 20 || resid < oldresid)) {
                        
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
        clean = resid < tol*wx*wy && ax >= 0. && ax <= wx && ay >= 0. && ay <= wy;

        if (resid > tol*wx*wy) {
#if 0
            // Try Lagrange multiplier method
            double lambda = 0.;
            double Phi = sqr(ax-x)+sqr(ay-y)+lambda*F(ax,ay);
            double Phix = 2.*(ax-x)+lambda*F(ax,ay);
            double Phiy = 2.*(ay-y)+lambda*F(ax,ay);
            double Phil = F(ax,ay);
            while (sqrt(sqr(Phix)+sqr(Phiy)+sqr(Phil)) > tol*wx*wy) {
                double ds = 1.;
                tax = ax-ds*Phix;
                tay = ay-ds*Phiy;
                double tal = lambda-ds*Phil;
                while (sqr(tax-x)+sqr(tay-y)+tal*F(tax,tay) > Phi) {
                    ds = ds/2.;
                    tax = ax-ds*Phix;
                    tay = ay-ds*Phiy;
                    tal = lambda-ds*Phil;
                }
                ax = tax;
                ay = tay;
                lambda = tal;
                Phi = sqr(ax-x)+sqr(ay-y)+lambda*F(ax,ay);
                Phix = 2.*(ax-x)+lambda*F(ax,ay);
                Phiy = 2.*(ay-y)+lambda*F(ax,ay);
                Phil = F(ax,ay);
            }
#else
#if 0
            // last ditch try bisection method
            double ax2;
            double ay2;
            double ax3, ay3;

            // if oscillating, it may be skipping over minimum
            // do one more iteration to get the other side
      
            p = F(ax,ay);
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            delta[0] = -p*px/(px*px+py*py);
            delta[1] = -p*py/(px*px+py*py);
            tax = ax+delta[0];
            tay = ay+delta[1];
                        
            p = F(tax,tay);
            px = Dx(tax,tay);
            py = Dy(tax,tay);
            double denom = px*px+py*py;
            delta[2] = py*(px*(tay-y)-py*(tax-x))/denom;
            delta[3] = -px*(px*(tay-y)-py*(tax-x))/denom;
                        
            if ((ax >= 0. && ax <= wx && ay >= 0. && ay <= wy && fabs(F(ax,ay))<tol)
                || (ax2 >= 0. && ax2 <= wx && ay2 >= 0. && ay2 <= wy && fabs(F(ax2,ay2))<tol)) {
                // compute cross product at each end
                double cp = Dx(ax,ay)*(ay-y)-Dy(ax,ay)*(ax-x);
                double cp2 = Dx(ax2,ay2)*(ay2-y)-Dy(ax2,ay2)*(ax2-x);
                if (cp*cp2 < 0.) {
                    // OK, we're oscillating, try bisection
                    double ax3 = (ax+ax2)/2.;
                    double ay3 = (ay+ay2)/2.;
                    while (fabs(F(ax3,ay3)) > tol*wx*wy 
                           && fabs(Dx(ax3,ay3)*(ay3-y)-Dy(ax3,ay3)*(ax3-x)) > tol*wx*wy) {
                        while (fabs(F(ax3,ay3)) > tol*wx*wy) {
                            p = F(ax3,ay3);
                            px = Dx(ax3,ay3);
                            py = Dy(ax3,ay3);
                            ax3 -= p*px/(px*px+py*py);
                            ay3 -= p*py/(px*px+py*py);
                        }
                        double cp3 = Dx(ax3,ay3)*(ay3-y)-Dy(ax3,ay3)*(ax3-x);
                        if (cp3*cp < 0.) {
                            ax2 = ax3;
                            ay2 = ay3;
                        } else {
                            ax = ax3;
                            ay = ay3;
                            cp = cp3;
                        }
                        ax3 = (ax+ax2)/2.;
                        ay3 = (ay+ay2)/2.;
                    }
                    ax = ax3;
                    ay = ay3;
                    clean = ax3 >= 0. && ax3 <= wx && ay3 >= 0. && ay3 <= wy;
                    if (clean) {
                        ax = ax3;
                        ay = ay3;
                    }
                }
            }
#else
	    return LocalDistByBoxes(x, y, ax, ay);
#endif
#endif
        }
        
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }

    double Bicubic::LocalDistNewton(const double x, const double y, 
                                    double& ax, double& ay, char& clean) const
    {
        ax = x;
        ay = y;
        double tol = 1.e-3;
        double p, px, py, pxx, pxy, pyy;
        double Jinv[2][2];
        double det;
        double delta[2];
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
                
        while (resid > tol*wx*wy && (itcount < 20 || resid < oldresid)) {
                        
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
                                                
            ax += delta[0];
            ay += delta[1];
            ++itcount;
            oldresid = resid;
            resid = sqrt(sqr(delta[0])+sqr(delta[1]));
        }
                
        //    clean = itcount < 20;
        clean = resid < tol*wx*wy && ax >= 0. && ax <= wx && ay >= 0. && ay <= wy;
        
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }

    class DistBox {
    public:
        double lx, ly;
        double dx, dy;
        double val[2][2];

        DistBox(void) : lx(0.), ly(0.), dx(0.), dy(0.) {}
        DistBox(const double x, const double y, const double delx, const double dely,
                const double v00, const double v10, const double v01, const double v11) 
            : lx(x), ly(y), dx(delx), dy(dely) 
            {val[0][0] = v00; val[1][0] = v10; val[0][1] = v01; val[1][1] = v11;}
    };
        
    double Bicubic::LocalDistByBoxes(const double x, const double y, 
                                     double& ax, double& ay) const
    {
        std::stack<DistBox> boxlist;
        int corner = (ax == 0 ? 0 : 1)+(ay == 0 ? 0 : 2);
        DistBox aBox(0., 0., wx, wy, F(0.,0.), F(wx,0.), F(0.,wy), F(wx,wy));
        DistBox bBox;
        boxlist.push(aBox);
        double dist = sqrt(wx*wx+wy*wy);
        double tol = 1.e-3*wx*wy;
        while (boxlist.size() > 0) {
            aBox = boxlist.top();
            boxlist.pop();
            if (sqrt(sqr(x-aBox.lx)+sqr(y-aBox.ly)) <= dist) {
                int score = (aBox.val[0][0]>0 ? 1 : 0)
                    +(aBox.val[0][1]>0 ? 2 : 0)
                    +(aBox.val[1][0]>0 ? 4 : 0)
                    +(aBox.val[1][1]>0 ? 8 : 0);
                if (score > 0 && score < 15) {
                    if (aBox.dx > tol || aBox.dy > tol) {
                        double p[3][3];
                        p[0][0] = aBox.val[0][0];
                        p[2][0] = aBox.val[1][0];
                        p[0][2] = aBox.val[0][1];
                        p[2][2] = aBox.val[1][1];
                        p[1][0] = F(aBox.lx+aBox.dx/2., aBox.ly);
                        p[0][1] = F(aBox.lx, aBox.ly+aBox.dy/2.);
                        p[1][1] = F(aBox.lx+aBox.dx/2., aBox.ly+aBox.dy/2.);
                        p[2][1] = F(aBox.lx+aBox.dx, aBox.ly+aBox.dy/2.);
                        p[1][2] = F(aBox.lx+aBox.dx/2., aBox.ly+aBox.dy);
                        switch (corner) {
                        case 0:
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[1][1], p[2][1], p[1][2], p[2][2]));
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[1][0], p[2][0], p[1][1], p[2][1]));
                            boxlist.push(DistBox(aBox.lx, aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[0][1], p[1][1], p[0][2], p[1][2]));
                            boxlist.push(DistBox(aBox.lx, aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[0][0], p[1][0], p[0][1], p[1][1]));
                            break;
                        case 1:
                            boxlist.push(DistBox(aBox.lx, aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[0][1], p[1][1], p[0][2], p[1][2]));
                            boxlist.push(DistBox(aBox.lx, aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[0][0], p[1][0], p[0][1], p[1][1]));
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[1][1], p[2][1], p[1][2], p[2][2]));
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[1][0], p[2][0], p[1][1], p[2][1]));
                            break;
                        case 2:
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[1][0], p[2][0], p[1][1], p[2][1]));
                            boxlist.push(DistBox(aBox.lx, aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[0][0], p[1][0], p[0][1], p[1][1]));
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[1][1], p[2][1], p[1][2], p[2][2]));
                            boxlist.push(DistBox(aBox.lx, aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[0][1], p[1][1], p[0][2], p[1][2]));
                            break;
                        case 3:
                            boxlist.push(DistBox(aBox.lx, aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[0][0], p[1][0], p[0][1], p[1][1]));
                            boxlist.push(DistBox(aBox.lx, aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[0][1], p[1][1], p[0][2], p[1][2]));
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly, aBox.dx/2., aBox.dy/2.,
                                                 p[1][0], p[2][0], p[1][1], p[2][1]));
                            boxlist.push(DistBox(aBox.lx+aBox.dx/2., aBox.ly+aBox.dy/2., aBox.dx/2., aBox.dy/2.,
                                                 p[1][1], p[2][1], p[1][2], p[2][2]));
                            break;
                        }
                    } else {
                        double tempdist = sqrt(sqr(x-aBox.lx-aBox.dx/2.)
                                               +sqr(y-aBox.ly-aBox.dy/2.));
                        if (tempdist < dist) {
                            dist = tempdist;
                            ax = aBox.lx+aBox.dx/2.;
                            ay = aBox.ly+aBox.dy/2.;
                        }
                    }
                }
            }
        }
        return dist;
    }

        
    double Bicubic::XCrossing(const double x, const double y, double& ay, char& clean) const
    {
        ay = y;
        double tol = 1.0e-3*wy;
        double delta = 1000.;
        int cnt = 0;
        while (fabs(delta) > tol && cnt < 20) {
            delta = -F(x,ay)/Dy(x,ay);
            ay += delta;
            ++cnt;
        }
        clean = fabs(delta) <= tol;
        return ay;
    }
        
    double Bicubic::YCrossing(const double x, const double y, double& ax, char& clean) const
    {
        ax = x;
        double tol = 1.0e-3*wx;
        double delta = 1000.;
        int cnt = 0;
        while (fabs(delta) > tol && cnt < 20) {
            delta = -F(ax,y)/Dx(ax,y);
            ax += delta;
            ++cnt;
        }
        clean = fabs(delta) <= tol;
        return ax;
    }

    void Bicubic::FollowToBdry(const int dir, const double x, const double y, double& ax, double& ay) const
    {
        ax = x; ay = y;
        double dt = 0.01;
        double tx[3], ty[3];
        double dp[2]; double p;
        double tdp[2]; double tp;
        double denom[2];
        while(ax > 0. && ax < wx && ay > 0. && ay < wy) {
            // take one step
            p = F(ax, ay);
            dp[0] = Dx(ax, ay);
            dp[1] = Dy(ax, ay);
            denom[0] = sqrt(dp[0]*dp[0]+dp[1]*dp[1]);
            double tol = 1;
            bool passed = false;
      
            while (!passed) {
                tx[0] = ax + dt*(dir*dp[1]/denom[0] - 10*p*dp[0]);
                ty[0] = ay + dt*(-dir*dp[0]/denom[0] - 10*p*dp[1]);
      
                // take two steps
                tx[1] = ax + dt/2.*(dir*dp[1]/denom[0] - 10*p*dp[0]);
                ty[1] = ay + dt/2.*(-dir*dp[0]/denom[0] - 10*p*dp[1]);
                tp = F(tx[1],ty[1]);
                tdp[0] = Dx(tx[1],ty[1]);
                tdp[1] = Dy(tx[1],ty[1]);
                denom[1] = sqrt(tdp[0]*tdp[0]+tdp[1]*tdp[1]);
                tx[2] = tx[1] + dt/2.*(dir*tdp[1]/denom[1] - 10*tp*tdp[0]);
                ty[2] = ty[1] + dt/2.*(-dir*tdp[0]/denom[1] - 10*tp*tdp[1]);
        
                double err = fabs(tx[0]-tx[2])+fabs(ty[0]-ty[2]);
                if (err > tol*dt)
                    dt = 0.9*dt*dt*tol/err;
                else {
                    passed = 1;
                    if (err < tol*dt/2)
                        dt = dt*2;
                }
            }

            ax = tx[0]; ay = ty[0];
        }
    }

    double Bicubic::Integral(void)
    {
        double ans = 0.;
        double wxp[4];
        double wyp[4];
        wxp[0] = wx;
        wyp[0] = wy;
        for (int i=1; i<4; ++i) {
            wxp[i] = wxp[i-1]*wx;
            wyp[i] = wyp[i-1]*wy;
        }
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                ans += a[i][j]/(i+1)/(j+1)*wxp[i]*wyp[j];

        return ans;
    }
}

